# ALAMO 1D Build and Support Guide

## Summary

In the current checkout, `./configure --dim=1 --comp=clang++` followed by
`make -j1` builds `bin/alamo-1d-clang++` as a **flame-only 1D executable**.
The 1D Makefile branch excludes the other integrator sources and creates only
the `alamo` launcher; `src/alamo.cc` accepts only `alamo.program=flame` in 1D.
The 2D/3D build paths remain unchanged.

This does not mean every ALAMO program or solid-mechanics model supports 1D.
Elastic mechanics, Sphere IC, and Ellipse IC are explicitly unavailable in the
1D flame build and abort if selected. A successful link is also not a validation
of a converted input or its physical equivalence.

The current 1D flame target has been built successfully. A two-step startup
smoke run also completed with thermal physics disabled and unused legacy inputs
allowed; see the test limitations below.

## What “1D support” needs to include

There are three distinct goals:

1. **Compile 1D executables.** All translation units needed by ALAMO must compile
   with `AMREX_SPACEDIM=1`.
2. **Run the desired physics in 1D.** Its integrator, material models,
   initial-condition builders, boundary conditions, and input parsing must all
   behave correctly in 1D.
3. **Validate and analyze 1D results.** Regression tests and any post-processing
   used for those results must accept 1D plotfiles.

Achieving only the first goal does not establish that a 1D pressure-oscillation
simulation is physically supported.

## Confirmed compile blockers

### 1. `Set::Vector` and `Set::Matrix` become the same type

[`Set/Base.H`](./src/Set/Base.H) defines `Set::Vector` as an Eigen matrix of shape
`AMREX_SPACEDIM x 1`, and `Set::Matrix` as one of shape
`AMREX_SPACEDIM x AMREX_SPACEDIM`. In 1D, both shapes are `1 x 1`, so both aliases
name the same C++ type. Overloads and explicit specializations that distinguish
them in 2D/3D then collide. The build output shows this in:

- [`String.H`](./src/Util/String.H), for `Join`
- [`Set.H`](./src/Set/Set.H), for `Field` specializations
- [`InputScraper.H`](./src/IO/InputScraper.H), for `InputArrayString`
- [`ParmParse.H`](./src/IO/ParmParse.H), for `queryarr` and related overloads

The fixes below describe a targeted approach: let the vector overload handle
one-value 1D representations in general utilities and parsers, while giving
vector and matrix fields a separate semantic role where their names and handling
must remain distinct.

#### What `String.H` and `Join` do

[`String.H`](./src/Util/String.H) is a collection of general-purpose string
helpers in the `Util::String` namespace. It provides operations such as replacing
text, wrapping lines, splitting strings, and joining a sequence of values into
one string. `Join` concatenates elements and inserts a separator between them;
it does not perform numerical calculations.

Some ALAMO call sites use `Join` to turn a vector or matrix default into the
space-separated text stored in an input parameter. Other uses join environment
or parameter values with underscores for output filenames, or join names for
diagnostic messages. For example, [`ParmParse.H`](./src/IO/ParmParse.H) uses it
when it inserts default array values into the parameter database, and
[`FileNameParse.cpp`](./src/IO/FileNameParse.cpp) uses it when replacing a
filename wildcard with an underscore-separated value.

For `Set::Vector` and `Set::Matrix`, a narrow 1D fix is to retain the vector
overload in every dimension and compile the matrix overload only in 2D/3D:

```cpp
AMREX_FORCE_INLINE
std::string Join(const Set::Vector & vec, char separator = ' ')
{
    // Existing vector implementation.
}

#if AMREX_SPACEDIM > 1
AMREX_FORCE_INLINE
std::string Join(const Set::Matrix & vec, char separator = ' ')
{
    // Existing matrix implementation.
}
#endif
```

In 1D, the matrix alias is accepted by the vector overload because it is exactly
the same type and contains one value. In 2D/3D the existing, distinct overloads
remain. The same guard pattern applies to `InputArrayString` in
[`InputScraper.H`](./src/IO/InputScraper.H): keep the vector implementation
unconditional and guard the matrix implementation with
`#if AMREX_SPACEDIM > 1`.

#### `ParmParse` array overloads

[`ParmParse.H`](./src/IO/ParmParse.H) has vector and matrix overloads for
`queryarr` and `queryarr_default`. Its matrix parser also contains compatibility
handling for 2x2 and 3x3 matrix inputs, so avoid replacing it with a generic
implementation that discards those behaviors.

Keep the vector overloads available in all dimensions and guard the matrix
overloads:

```cpp
int queryarr(std::string name, Set::Vector & value, Unit unit = Unit::Less(),
             const std::source_location &location =
                 std::source_location::current())
{
    // Existing vector parser.
}

#if AMREX_SPACEDIM > 1
int queryarr(std::string name, Set::Matrix & value, Unit unit = Unit::Less(),
             const std::source_location &location =
                 std::source_location::current())
{
    // Existing matrix parser, including 2D/3D compatibility handling.
}
#endif
```

Apply the same guard to the `queryarr_default` overload that takes
`Set::Matrix` and a matrix default. In 1D, calls using `Set::Matrix` resolve to
the vector overload, so both aliases follow the same one-value path. Review the
vector parser's length policy: it currently accepts at least
`AMREX_SPACEDIM` values and copies the first dimensions. Decide whether extra
values should be rejected in 1D, and test that behavior rather than
unintentionally changing it while fixing the overload collision.

#### `Field` specializations need a semantic role

[`Set.H`](./src/Set/Set.H) specializes `Field` operations for vector and matrix
values. In 1D, `Field<Set::Vector>` and `Field<Set::Matrix>` are the same C++
specialization. The `Copy` and `NComp` operations both represent one stored
component, but `Name` must preserve the difference between a vector component
such as `_x` and a matrix component such as `_xx`. Simply removing the duplicate
matrix specialization could compile while giving matrix output the wrong name.

A robust approach is to add a role parameter to `Field`, keeping the stored
Eigen value type as `T`. Make these edits in [`Set.H`](./src/Set/Set.H):

1. Immediately before the existing `template <class T> class Field`
   declaration, add the two role tags. Replace the existing template line; do
   **not** leave both template declarations one after another. The result should
   look like this:

   ```cpp
   struct VectorFieldRole {};
   struct MatrixFieldRole {};

   template <class T, class Role = void>
   class Field : public amrex::Vector<
       std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>
   {
       // Keep the existing Field implementation here.
   };
   ```

   In context, the beginning of the class section should be:

   ```cpp
   using HC = Hypercube;

   struct VectorFieldRole {};
   struct MatrixFieldRole {};

   template <class T, class Role = void>
   class Field : public amrex::Vector<
       std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>
   {
       // Existing constructors and methods.
   };
   ```

   If you already inserted a standalone `template <class T, class Role = void>`
   above the old `template <class T>`, delete the standalone line and replace
   the old line with the role-parameter version. A `template` line must directly
   precede the class or function it parameterizes.

2. Immediately after the closing `};` of `Field` and before the existing
   `template <typename T> using Patch = ...` alias, add:

   ```cpp
   using VectorField = Field<Set::Vector, VectorFieldRole>;
   using MatrixField = Field<Set::Matrix, MatrixFieldRole>;
   ```

   These aliases name field types only; do not replace `Set::Vector` or
   `Set::Matrix` in Eigen-valued calculations.

3. Keep the existing method bodies. Change the type in every vector
   specialization from `Field<Set::Vector>` to
   `Field<Set::Vector, VectorFieldRole>`. In the current file, that means the
   `Copy`, `CopyFrom`, `Add`, `AddFrom`, `NComp`, and `Name` specializations.
   For example, change:

   ```cpp
   void Field<Set::Vector>::Copy(
       int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
   ```

   to:

   ```cpp
   void Field<Set::Vector, VectorFieldRole>::Copy(
       int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
   ```

   Apply that same type change to `CopyFrom`, `Add`, and `AddFrom`. For the two
   one-line vector specializations, keep the bodies and change the signatures
   exactly as follows:

   ```cpp
   template<> ALAMO_SINGLE_DEFINITION
   int Field<Set::Vector, VectorFieldRole>::NComp() const
   {
       return AMREX_SPACEDIM;
   }

   template<> ALAMO_SINGLE_DEFINITION
   std::string Field<Set::Vector, VectorFieldRole>::Name(int i) const
   {
       // Keep the existing vector Name body, which returns name + "_x" in 1D.
   }
   ```

   In the current one-line `NComp` definition, this is specifically the
   replacement:

   ```cpp
   // Before
   template<> ALAMO_SINGLE_DEFINITION int
   Field<Set::Vector>::NComp() const {return AMREX_SPACEDIM;}

   // After
   template<> ALAMO_SINGLE_DEFINITION int
   Field<Set::Vector, VectorFieldRole>::NComp() const
   {
       return AMREX_SPACEDIM;
   }
   ```

   For the vector `Name` definition, the only required change is the class
   template argument; leave its component-name logic unchanged:

   ```cpp
   // Before
   template<> ALAMO_SINGLE_DEFINITION std::string
   Field<Set::Vector>::Name(int i) const

   // After
   template<> ALAMO_SINGLE_DEFINITION std::string
   Field<Set::Vector, VectorFieldRole>::Name(int i) const
   ```

   For the matrix specializations (`Copy`, `NComp`, and `Name`), make the
   corresponding substitutions:

   ```cpp
   // Before
   Field<Set::Matrix>::Copy(...)
   Field<Set::Matrix>::NComp()
   Field<Set::Matrix>::Name(...)

   // After
   Field<Set::Matrix, MatrixFieldRole>::Copy(...)
   Field<Set::Matrix, MatrixFieldRole>::NComp()
   Field<Set::Matrix, MatrixFieldRole>::Name(...)
   ```

   The matrix `NComp` body remains `return AMREX_SPACEDIM *
   AMREX_SPACEDIM;`. The matrix `Name` body remains unchanged and returns
   `name + "_xx"` for component zero in 1D. `NComp` is one for both fields in
   1D; the role is needed so the separate `Name` specializations still compile
   and preserve the distinct labels.

4. Update field declarations and function parameters elsewhere in `src/`.
   Replace the exact field-type spellings:

   ```cpp
   Set::Field<Set::Vector>  // becomes Set::VectorField
   Set::Field<Set::Matrix>  // becomes Set::MatrixField
   ```

   For example:

   ```cpp
   Set::Field<Set::Vector> disp_mf;
   Set::Field<Set::Matrix> stress_mf;
   ```

   becomes:

   ```cpp
   Set::VectorField disp_mf;
   Set::MatrixField stress_mf;
   ```

   Make the same substitutions in parameters and references, including
   `Set::Field<Set::Vector>&`. Do not change `Set::Field<Set::Scalar>` or
   fields of other types. The current vector/matrix field uses include
   `Integrator/Base/Mechanics.H`, `Integrator/LowMach.H`,
   `Solver/Nonlocal/Newton.H`, `IC/Constant.H`, `IC/Expression.H`, `IC/Trig.H`,
   and `BC/Operator/Elastic/Elastic.H`.

   You can do these edits across the source folder in VS Code:

   1. Open **Search and Replace** with `Ctrl+Shift+H`.
   2. Turn **Use Regular Expression** off so the angle brackets are treated
      literally.
   3. In **files to include**, enter `src/**`.
   4. Search for `Set::Field<Set::Vector>` and replace with
      `Set::VectorField`. Review the matches, then choose **Replace All**.
   5. Search for `Set::Field<Set::Matrix>` and replace with
      `Set::MatrixField`. Review the matches, then choose **Replace All**.

   These exact searches change field declarations and parameters, not ordinary
   Eigen value types such as `Set::Vector`. They also do not match the
   specializations inside `Set.H`, which start with `Field<...>` and must be
   updated separately as described in step 3.

   Then search for remaining old field uses with:

   ```bash
   rg 'Set::Field<Set::(Vector|Matrix)>' src
   ```

   Convert matches that are field declarations, parameters, or references.
   The class specializations inside `Set.H` are handled separately in step 3.
   Do not blindly replace every `Set::Vector` or `Set::Matrix`; most of those
   are value types and should stay unchanged.

The role parameter makes the vector and matrix **field types** distinct even
when their stored Eigen **value type** is the same 1x1 matrix. This avoids
wrapping Eigen and preserves the component labels in output.

The role tag is a compile-time type distinction; it adds no data member and
does not change the Eigen value, component count, copy arithmetic, or mesh
storage. In 2D/3D the existing vector and matrix value types are already
distinct, so the role should not change numerical results when the existing
specialization bodies and call sites are preserved. However, the C++ field
types (and their mangled symbols) do change for vector and matrix fields in
every dimension. Update all signatures and declarations that use those fields
and rebuild every executable; do not mix objects or libraries built before and
after the template change. Run the existing 2D tests and compare representative
outputs to confirm behavior was preserved. Apply equivalent checks to 3D if
that build is part of the supported workflow.

### 2. Scalar boundary-condition classes omit 1D faces

[`BC/Constant.H`](./src/BC/Constant.H) and
[`BC/Expression.H`](./src/BC/Expression.H) define face enums, face counts, and
face-name tables for 2D and 3D, but not 1D. Their later code still uses those
definitions to size arrays, parse boundary-condition keys, and access faces.
This causes the reported undeclared `Face`, `m_nfaces`, and `facestr` errors.
Some paths also refer directly to `YLO`/`YHI`.

Both classes need a 1D enum with two active faces and an `INT` sentinel. Add the
following dimension branch to each class's existing enum:

```cpp
#if AMREX_SPACEDIM == 1
enum Face { XLO, XHI, INT };
#elif AMREX_SPACEDIM == 2
// Existing 2D enum.
#elif AMREX_SPACEDIM == 3
// Existing 3D enum.
#endif
```

Their face arrays also need a compile-time size of two in 1D:

```cpp
#if AMREX_SPACEDIM == 1
static constexpr int m_nfaces = 2;
#elif AMREX_SPACEDIM == 2
static constexpr int m_nfaces = 4;
#elif AMREX_SPACEDIM == 3
static constexpr int m_nfaces = 6;
#endif
```

The first 1D build also exposed dimension assumptions in
[`BC/Constant.cpp`](./src/BC/Constant.cpp), the implementation behind
`Constant.H`:

- The legacy constructor initializes y-face types and values and accesses
  `_bc_lo_2`/`_bc_hi_2`. Those parameter names are removed by `AMREX_D_DECL`
  in 1D. Guard the y-face initialization and those parameter accesses with
  `#if AMREX_SPACEDIM > 1`.
- `FillBoundary` loads y-face type/value locals and contains ylo/yhi branches
  in its `ParallelFor` lambda. Guard both the locals and the y-face branches
  with `#if AMREX_SPACEDIM > 1`. The z-face branch already has a 3D guard.
- The corner-averaging kernel after the boundary-fill kernel assumes 2D
  corners, as its comment notes. Guard it for dimensions greater than one or
  replace it with a dimension-general implementation before using it in 1D.

In [`Expression.H`](./src/BC/Expression.H), put the name table inside the
`BC::Expression` class, immediately after `static constexpr const char* name =
"expression";` and before the `enum Face` declaration. Replace the existing
dimension-conditional `facestr` block there; do not put this table inside
`Parse()` or next to its face loop. The 1D branch should be:

```cpp
#if AMREX_SPACEDIM == 1
static constexpr const char * const facestr[] = {"xlo", "xhi"};
#elif AMREX_SPACEDIM == 2
// Existing 2D table.
#elif AMREX_SPACEDIM == 3
// Existing 3D table.
#endif
```

These definitions resolve the missing names and let `Expression::Parse`'s
`for (face = 0; face != Face::INT; ++face)` loop visit `xlo` and `xhi`. They do
not make the explicit y-face indexing in either parser safe. In `Constant.H`,
the 1D `Face` enum is not enough on its own. Also add
`static constexpr int m_nfaces = 2` to the dimension-conditional block where
the 2D and 3D face counts are declared. In `Constant::Parse`, replace the
explicit x/y/z face-clearing statements with a loop over `m_nfaces` (the exact
code is in the targeted edit steps below). That avoids referring to inactive
`YLO`/`YHI` enum values in 1D. Guard the separate y-face type and value parsing
blocks with `#if AMREX_SPACEDIM > 1`; those blocks index `Face::YLO`/`Face::YHI`.
In the 1D `#else`, ignore `type.ylo`, `type.yhi`, `val.ylo`, and `val.yhi`.
The 3D z-value blocks are already dimension-guarded; the `type.zlo`/`type.zhi`
keys are queried by the existing parser.

In `Expression::Parse`, find the `// TYPES` comment and the two consecutive
blocks that start at:

```cpp
pp.queryarr_default("type.ylo", str, {"dirichlet"});
```

and

```cpp
pp.queryarr_default("type.yhi", str, {"dirichlet"});
```

They come immediately after the `type.xhi` block and before the
`pp.queryarr_default("type.zlo", ...)` block. Put `#if AMREX_SPACEDIM > 1`
immediately before the existing `str.clear()` that precedes the `type.ylo`
call. Keep both existing y-face blocks and the `type.zlo` query inside that
branch. After the `type.zlo` query, add the `#else` and the ignore calls below;
close with `#endif` before the existing `#if AMREX_SPACEDIM == 3` block that
processes the zlo face. The resulting structure is:

```cpp
#if AMREX_SPACEDIM > 1
pp.queryarr_default("type.ylo", str, {"dirichlet"});
// Keep the existing type.ylo parsing code.
pp.queryarr_default("type.yhi", str, {"dirichlet"});
// Keep the existing type.yhi parsing code.
pp.queryarr_default("type.zlo", str, {"dirichlet"});
#else
pp.ignore("type.ylo");
pp.ignore("type.yhi");
pp.ignore("type.zlo");
#endif
```

Here, the comments stand for the existing validation, indexing, and resize
statements—do not delete those statements. Put this guard in the type-parsing
section, not around the later `for (face = 0; face != Face::INT; ++face)` value
loop. In 1D that loop already visits only `xlo` and `xhi`; the `pp.ignore`
calls mark supplied inactive face keys as consumed so ALAMO does not later
report them as unused. The `type.zhi` query is outside that block and remains
consumed by the existing code; do not both ignore and query `type.zhi`.

If the input format supports combined face-key variants, add them to the
lower-dimensional ignore policy too. Finally, verify `GetBCTypes()`'s high-side
row while editing: the current code appears to select low-face entries for both
rows. If the rows represent low and high sides, the second row should use the
`XHI`/`YHI`/`ZHI` entries. Test both sides rather than relying on compilation
alone.

### 3. Elastic boundary-condition operators also omit 1D

[`BC/Operator/Elastic/Elastic.H`](./src/BC/Operator/Elastic/Elastic.H) defines its
face enum only for 2D/3D. Its constant and expression implementations consume that
enum, so this is another compile blocker for a complete ALAMO build. It may not
be exercised by the flame case, but it must be handled if those translation
units are part of the build.

### 4. The `Matrix4::Major` 1D path is incomplete

[`Matrix4_Major.H`](./src/Set/Matrix4_Major.H) lacks a complete 1D implementation
for operations instantiated by ALAMO's tests. The test executable instantiates
these matrix operations in [`test.cc`](./src/test.cc). Implement and test the
1D behavior, or explicitly define which operations are unsupported and ensure
the build/test configuration does not instantiate them accidentally.

### Full 1D support remains out of scope for the flame-only build

The IC interface carries `Set::VectorField` through a `VectorIC` base, so
vector `Add` overrides remain available in 1D. `Sphere` and `Ellipse` compile
but abort with an explicit unsupported-feature message if selected. Elastic
`Constant`, `Expression`, and `TensionTest` have the same policy.

A previous build that compiled all integrators exposed broader solid-mechanics
assumptions. Those paths are excluded from the current flame-only 1D build, but
would need work before claiming full ALAMO 1D support:

- [`Matrix4_Major.H`](./src/Set/Matrix4_Major.H) lacks 1D indexing,
  subtraction, and norm operations used by solid-model checks in
  [`Solid.H`](./src/Model/Solid/Solid.H).
- The earlier vector/matrix overload collisions in
  [`Numeric/Stencil.H`](./src/Numeric/Stencil.H) are now guarded for active
  dimensions greater than one. Keep checking remaining call sites as LowMach
  is enabled; this does not validate the full stencil library in 1D.
- [`Mechanics.H`](./src/Integrator/Base/Mechanics.H) expects a 3x3 `Matrix3`
  from a gradient that returns the 1x1 `Set::Matrix` in 1D.
- [`Model/Solid/Linear/Cubic.H`](./src/Model/Solid/Linear/Cubic.H) requests
  2D/3D `Matrix4::Cubic` operations and fourth-order tensor products not
  implemented for dimension 1.

The flame-only configuration avoids those unsupported modules; it does not
implement or validate them for 1D.

### LowMach/Rocfire 1D port (in progress)

The default 1D build remains the flame-only `alamo-1d-clang++`. An opt-in
`make lowmach-1d -j1` target now builds `bin/lowmach-1d-clang++` without
changing that default. The first supported scope is LowMach with rigid
condensed phases; deformable-solid mechanics are explicitly rejected in 1D.
This keeps the Beckstead AP use case in scope without pretending the
3D-based NeoHookean model is a 1D implementation.

The port also carries field-role tags through `Integrator::AddField` and its
AMR wrapper, makes LowMach diagnostics and timestep speed calculations use
active dimensions, and excludes the finite-solid stress path from 1D. The
1D nodal AMR path uses the NodeBilinear interpolator, and its mapper has a
persistent lifetime for AMReX's fill-patch cache.

Validation so far:

- `make lowmach-1d -j1` builds the new executable; the default
  `make -j1` still builds `alamo-1d-clang++`.
- A focused compile of `LowMach.cpp` against the 2D AMReX headers succeeds.
  This is not a full 2D link or regression run.
- [`input.lm.rocfire_1d_smoke`](./input.lm.rocfire_1d_smoke) completed a
  five-step, strict-parser 1D LowMach/Rocfire startup. It uses 64 cells on
  `[-0.15, 0.05] mm`, rigid AP, a short heat pulse, and constant 3 MPa
  pressure; all sampled fields were finite.
- The same setup, extended by command-line overrides to about 1.5 ms with
  dynamic timestepping and implicit viscosity/thermal/species diffusion,
  completed 1,510 steps without NaNs. AP gas and Rocfire intermediates became
  nonzero; the final maximum temperature was about 1,107 K. This confirms
  execution of the 1D chemistry/phase-change path, not a converged burn rate.
  Plotfiles are under
  `results/lowmach_1d_rocfire_smoke/chem_1D_20261109_145722`.
  The run used:

  ```bash
  ./bin/lowmach-1d-clang++ input.lm.rocfire_1d_smoke \
    timestep=1.0e-7_s stop_time=1.5e-3_s max_step=15000 \
    dynamictimestep.on=1 dynamictimestep.min=1.0e-12 \
    dynamictimestep.max=1.0e-6 viscosity=2 thermal_diffusion=2 \
    species_diffusion=2 amr.plot_int=-1 amr.plot_dt=5.0e-5_s \
    plot_file=results/lowmach_1d_rocfire_smoke/chem_%DD_%Y%m%d_%H%M%S
  ```
- A five-step AMR smoke with `amr.max_level=1` created both levels and had no
  nonfinite temperature, pressure, velocity, AP-gas, or AP-solid values.
  Extended diagnostics also completed for five steps with every plotfile field
  finite. Outputs are under
  `results/lowmach_1d_rocfire_smoke/amr_1D_20261109_145819` and
  `results/lowmach_1d_rocfire_smoke/diagnostics_1D_20261109_150344`.

- [`input.lm.beckstead_1d_130hz_smoke`](./input.lm.beckstead_1d_130hz_smoke)
  adapts the 2D Beckstead pilot: 500 psia, 130 Hz, 2% pressure drive, the
  2-mm AP / 1-mm gas domain, and the calibrated pure-AP Rocfire/phase-change
  parameters. It maps the axial y direction to x, uses 240 base cells and
  AMR through level 2 (3.125 um finest spacing), and runs to 10 ms. The
  pressure drive starts at 1.75 ms, so this covers only about 1.07 driven
  cycles.
- A matched zero-drive case also completed to 10 ms. The eta=0.5 front moved
  3.863 um in the driven case versus 3.854 um in the control. Their final
  difference, 9.3 nm, is far below the 3.125 um finest cell; this run does not
  resolve a pressure-response signal. It is only evidence that the 1D
  Beckstead input parses and runs with chemistry, heat transfer, and pressure
  forcing enabled. Do not use it to fit a transfer function.
- The 10-ms forced/control plot is
  [`analysis_10ms_smoke.png`](./results/beckstead_1d_lowmach_130hz/analysis_10ms_smoke.png).
  Its nominal post-drive-window fitted rates are about 0.339 and 0.338 mm/s;
  that fit remains an ignition transient. The plot annotates the front
  difference against the finest-cell scale.
- VisIt’s installed readers do not include a 1D BoxLib reader, so its
  “Boxlib2D/Boxlib3D/Silo” error is expected for these valid 1D plotfiles. The
  new [`export_1d_plotfile.py`](./scripts/export_1d_plotfile.py) uses yt to
  export a leaf-cell `final_profile.csv` and an `front_history.csv` from a
  plot directory:

  ```bash
  python3 scripts/export_1d_plotfile.py \
    results/beckstead_1d_lowmach_130hz/output_1D_20261111_141507
  ```

- The forced/control outputs and CSV exports are under
  `results/beckstead_1d_lowmach_130hz/`. The 10-ms Beckstead smoke is
  under-resolved for response analysis: its maximum driven/control front
  difference is about 9.3 nm, versus a 0.78125-um finest cell. It is useful
  only as a setup check.
- A level-4 driven attempt with a 100-ms startup trim was launched to the
  existing 0.43462-s pilot stop time. It advanced to about 0.39402 s before
  the pressure-projection MLMG solver failed as the timestep collapsed. The
  partial history shows the front reaching x=-1.994 mm at 0.3938 s, nearly
  exhausting the 2-mm AP slab; AP-gas density and temperature then spike.
  Treat this as a domain-depth failure and discard the late-time data, not as
  a pressure-response result. The partial plot is
  [`analysis_partial.png`](./results/beckstead_1d_lowmach_130hz/trim100_driven_1D_20261111_142834/analysis_partial.png).
- A deeper matched pair completed successfully with the domain extended to
  x=[-3, 1] mm, 320 base cells, AMR level 4, and the same 0.43462-s stop time.
  Both cases completed 86,925 steps and produced 2,174 snapshots. At the end,
  the front was near x=-2.193 mm, leaving about 0.807 mm before the lower
  boundary. The 100-ms-trim histories span 0.33462 s (about 43.5 cycles at
  130 Hz, since the input drive continues through the stop time).
- A separate matched-geometry 1D pair also completed at the 2D pilot's
  2-mm-AP/1-mm-gas depth, 240 axial cells, level 4, and 0.28462-s stop time.
  The common 100-280 ms interval has nearly identical axial resolution and
  output cadence in the 1D and 2D data. A fixed-frequency fit to the
  driven-minus-control front history gives a 130-Hz burn-rate amplitude of
  0.0949 mm/s (1.91% of mean) in 2D and 0.1098 mm/s (2.07%) in 1D; fitted
  phases differ by about 2.5 degrees. Thus the 1D run does reproduce the
  oscillatory response seen in 2D. The exact-geometry comparison is plotted
  in
  [`compare_1d_2d_130hz_response_2mm.png`](./results/beckstead_1d_lowmach_130hz/compare_1d_2d_130hz_response_2mm.png)
  with fit values in
  [`compare_1d_2d_130hz_metrics_2mm.csv`](./results/beckstead_1d_lowmach_130hz/compare_1d_2d_130hz_metrics_2mm.csv).
  The companion
  [`compare_1d_2d_130hz_raw_response_and_means_2mm.png`](./results/beckstead_1d_lowmach_130hz/compare_1d_2d_130hz_raw_response_and_means_2mm.png)
  plots centered finite-difference responses without polynomial or harmonic
  smoothing, with the window-mean rates in a second panel.
- The absolute mean burn rate remains about 6.4% higher in 1D: 5.2933 mm/s
  versus 4.9770 mm/s in 2D. The same offset appears in the zero-drive controls,
  and the matched 2-mm 1D run agrees with the 3-mm 1D run, so neither the
  pressure forcing nor the deeper 3-mm boundary explains it. In 2D, the front
  spread across the periodic width is under 0.08 um, so front corrugation also
  does not explain the offset. The inputs have matching scalar physics and
  pressure forcing; this baseline-rate difference remains a dimensional
  implementation/build or measurement question.
- The earlier 1D rate plots used a 10-ms moving fit, while the 2D plot uses
  centered finite differences of front positions. A 10-ms window spans 1.3
  cycles at 130 Hz and attenuates the matched 2-mm 1D rate oscillation to
  about 6% of its centered-difference amplitude. Use centered differences or
  a phase-locked fit when visually comparing oscillations. The raw comparison
  plot differentiates the driven-minus-control front positions directly; it
  does not detrend or fit the response. Both data sets have more than 40
  samples per cycle.
- The front displacement difference grows from under 1 um at 20 ms to about
  10 um by 30 ms, then continues to grow at the faster 1D mean rate. This
  early drift is also present in the zero-drive comparison; the pressure
  oscillation is not its cause.
- A 3-mm level-5 1D pair provides a spatial-refinement check: relative to
  level 4, its 200-ms-trim mean rate changes by about 0.1%, the
  control-corrected response magnitude by about 0.3%, and phase by about
  0.12 degrees. The level-5 response diagnostics are in
  [`frequency_response.csv`](./results/beckstead_1d_lowmach_130hz/mesh5_3mm_30cycles_response_200ms/frequency_response.csv).
- The 1D oscillatory response is promising, but the mean-rate discrepancy,
  archived 2D solver-build provenance, and temporal convergence remain open
  validation items. Do not treat the current result as a validated transfer
  function. See the expanded [1D status and validation summary](./1D_STATUS_SUMMARY.md).
  The 3-mm pair comparison is
  [`trim100_3mm_30cycles_comparison.png`](./results/beckstead_1d_lowmach_130hz/trim100_3mm_30cycles_comparison.png);
  individual history and final-profile plots are in each run directory.
- For a single-case overview, run
  [`plot_1d_history.py`](./scripts/plot_1d_history.py) on an exported
  `front_history.csv`. For paired runs, use
  [`analyze_1d_regression.py`](./scripts/analyze_1d_regression.py); its plot
  includes front positions, local fitted regression rates, drive/control
  front difference, and pressure range.
- Do not use the single-grid response as a validated transfer function yet.
  First establish spatial and temporal convergence and reconcile the 6.4%
  baseline-rate difference with an otherwise matched 2D case.
- To visualize the final spatial state, use
  [`plot_1d_profile.py`](./scripts/plot_1d_profile.py) on the exported
  `final_profile.csv`; it plots temperature, pressure, rigid eta, AP solid,
  and AP gas through the interface.

## Targeted edit sequence for the current compiler errors

Use the following source anchors rather than editing every repeated compiler
diagnostic. These steps correspond to the current source layout; nearby comments
or function signatures are more reliable anchors than line numbers if the files
change.

### A. Guard matrix-only overloads in utility and input code

1. In [`String.H`](./src/Util/String.H), find the adjacent functions beginning
   with `Join(const Set::Vector & vec, ...)` and
   `Join(const Set::Matrix & vec, ...)`. Leave the vector function unchanged.
   Insert `#if AMREX_SPACEDIM > 1` immediately before the matrix function's
   `AMREX_FORCE_INLINE`, and `#endif` after that function's closing brace:

   ```cpp
   #if AMREX_SPACEDIM > 1
   AMREX_FORCE_INLINE
   std::string Join(const Set::Matrix & vec, char separator = ' ')
   {
       // Keep the existing matrix body.
   }
   #endif
   ```

   This removes the duplicate overload only in 1D. A 1D matrix call then binds to
   the vector overload and serializes its one element.

2. In [`InputScraper.H`](./src/IO/InputScraper.H), find
   `InputArrayString(const Set::Vector &value)` followed by
   `InputArrayString(const Set::Matrix &value)`. Leave the vector overload in
   place; wrap only the matrix overload, from its declaration through its
   closing brace, in the same `#if AMREX_SPACEDIM > 1` guard. The 1D vector
   serializer handles the equivalent one-element matrix.

3. In [`ParmParse.H`](./src/IO/ParmParse.H), guard the matrix-specific
   `queryarr` overload. Its signature starts
   `queryarr(std::string name, Set::Matrix & value, Unit unit ...)`; the function
   ends immediately before the `template<typename T> queryarr_required(...)`
   overload. Place `#if AMREX_SPACEDIM > 1` before the matrix signature and
   `#endif` after its closing brace.

   Write the preprocessor line exactly as `#if AMREX_SPACEDIM > 1`; do not
   append the explanatory word `guard`.

   Do the same for the matrix `queryarr_default` overload: start at
   `queryarr_default(std::string name, Set::Matrix & value, Set::Matrix def, ...)`
   and end before the following `queryarr_default` overload for
   `std::vector<double>`. Do not guard the vector overloads. In 1D, matrix calls
   will then resolve to the vector versions; in 2D/3D the existing matrix parser,
   including its 2x2/3x3 compatibility behavior, remains active.

### B. Keep vector fields and matrix fields semantically distinct

In [`Set.H`](./src/Set/Set.H), add `VectorFieldRole` and `MatrixFieldRole`
immediately before the `Field` class template. Change the class declaration from
`template <class T>` to `template <class T, class Role = void>`. After the class
definition, add the `VectorField` and `MatrixField` aliases shown in the
preceding section.

Next, update each explicit vector specialization (`Copy`, `CopyFrom`, `Add`,
`AddFrom`, `NComp`, and `Name`) to specialize
`Field<Set::Vector, VectorFieldRole>`. Update the matrix specializations
(`Copy`, `NComp`, and `Name`) to specialize
`Field<Set::Matrix, MatrixFieldRole>`. The role makes these distinct template
specializations even when the value type is the same Eigen 1x1 type.

Finally, update field declarations and function parameters throughout `src/`:
replace vector field types with `Set::VectorField` and matrix field types with
`Set::MatrixField`. Search before and after the migration:

```bash
rg 'Field<Set::(Vector|Matrix)>' src
```

Convert every use that refers to a vector or matrix field; leave fields of other
types, such as `Field<Set::Scalar>`, unchanged. This migration is needed to keep
the `_x` and `_xx` output names distinct. Merely wrapping the matrix member
specializations in `#if AMREX_SPACEDIM > 1` is not a complete fix: in 1D it
would make matrix fields use the vector field's name specialization.

### C. Add 1D faces and make parsing dimension-aware

In [`Constant.H`](./src/BC/Constant.H), add the 1D `Face` branch before its
existing 2D branch, and add the `m_nfaces = 2` branch before the existing
`m_nfaces = 4` branch. In [`Expression.H`](./src/BC/Expression.H), do the same
for its `Face` enum and `m_nfaces`, and add the 1D `facestr` table before its
existing 2D table. The required definitions are shown in the scalar boundary
condition section above.

Then update the parse code so it accesses only active faces:

- In `Constant::Parse`, replace the current x/y/z face-clearing statements with:

  ```cpp
  for (int face = 0; face < m_nfaces; ++face)
  {
      value.m_bc_type[face].clear();
      value.m_bc_val[face].clear();
  }
  ```

  This clears exactly the active face slots: two in 1D, four in 2D, six in 3D.
  Then guard the `type.ylo`/`type.yhi` block: put
  `#if AMREX_SPACEDIM > 1` before the existing `str.clear()` immediately before
  `pp.queryarr_default("type.ylo", ...)`. Keep both existing y type blocks and
  the `type.zlo` query inside it. Add an `#else` with
  `pp.ignore("type.ylo");`, `pp.ignore("type.yhi");`, and
  `pp.ignore("type.zlo");`. Close with `#endif` before the existing
  `#if AMREX_SPACEDIM == 3` block that processes the zlo face. The current
  misplaced type guard has a closing `#endif` at the end of `Parse()`; remove
  that old closing directive when you add this correctly placed one. Do not
  ignore `type.zhi`; it is queried later outside the guard.

  Do the same for the values: put `#if AMREX_SPACEDIM > 1` immediately before
  `value.m_bc_val[Face::YLO].clear()`, and include the existing ylo and yhi
  value handling through the last `else Util::Exception(...)` for `val.yhi`.
  Keep the `val.clear()` between the ylo and yhi blocks inside this guard too.
  Do not close the guard immediately after `value.m_bc_val[Face::YLO].clear()`;
  the `Face::YLO`/`Face::YHI` indexing is in the parsing code that follows.
  In its `#else`, ignore `val.ylo` and `val.yhi`; close with `#endif` only after
  the final error-handling line for `val.yhi`, before the existing `val.clear()`
  and `val.zlo` query. Keep xlo/xhi parsing active in all dimensions. Do not put
  the `val.zlo`/`val.zhi` queries inside this guard; their current code consumes
  those keys and only indexes z faces in the 3D-guarded blocks. This is a second,
  independent `#if`/`#else`/`#endif`; do not leave the type guard open across
  the value parser.
- In `Expression::Parse`, guard the y-face type handling beginning at
  `pp.queryarr_default("type.ylo", ...)` and ending after the yhi block, before
  `type.zlo`. In the 1D `#else`, ignore `type.ylo`, `type.yhi`, and inactive
  `val.*` keys. Its `for (face = 0; face != Face::INT; ++face)` value loop will
  automatically process only xlo/xhi once the 1D enum, face count, and name
  table are defined.
- For inactive y/z input keys, explicitly choose whether to reject or ignore
  them. If keeping lower-dimensional inputs backward-compatible, mark the
  inactive `type.*` and `val.*` keys ignored using `pp.ignore(...)`, following
  the existing 2D handling for unused z-face keys in `Expression::Parse`.
  Include any combined face-key variants that the input parser recognizes.

While editing each `GetBCTypes()` method, inspect its two returned sides. The
current code appears to select low-face entries for both rows. If these rows
represent low and high sides, use `XLO`/`YLO`/`ZLO` in the first row and
`XHI`/`YHI`/`ZHI` in the second. Add a test that would fail if either side is
incorrect.

### D. Make the `Constant.cpp` implementation dimension-aware

The current compiler errors in `src/BC/Constant.cpp` are from y-face
assumptions in the legacy constructor and `FillBoundary`:

- Guard the constructor's `Face::YLO`/`Face::YHI` initialization and its
  `_bc_lo_2`/`_bc_hi_2` references. Those parameter names are not declared in
  1D because `AMREX_D_DECL` removes them.
- In `FillBoundary`, guard the `bc_type_ylo`/`bc_type_yhi` and
  `bc_val_ylo`/`bc_val_yhi` locals, then guard the two y-face `else if`
  branches in the `ParallelFor` lambda with `#if AMREX_SPACEDIM > 1`.
- Guard the 2D corner-averaging `ParallelFor` after the boundary-fill kernel,
  or generalize it; its current logic assumes y-face corners.

After these changes, search `Constant.cpp` for `Face::YLO`, `Face::YHI`,
`Orientation::ylo`, `Orientation::yhi`, `_bc_lo_2`, and `_bc_hi_2`; each use
must be within a dimension guard or absent from the 1D preprocessing path.
`GetBCRec`, `IsPeriodic`, and the periodicity methods use `AMREX_D_DECL` for
dimension-dependent arguments, so inspect their 1D macro expansion before
adding extra guards.

The next serial build exposed the same implementation-level assumption in
[`BC/Expression.cpp`](./src/BC/Expression.cpp). In `Expression::FillBoundary`,
guard the y-face type/function locals and ylo/yhi branches with
`#if AMREX_SPACEDIM > 1`. Guard the corner-averaging kernel after the boundary
kernel as well; it assumes 2D corners. `GetBCRec`, `IsPeriodic`, and the
periodicity methods use `AMREX_D_DECL` for their dimension-dependent arguments.

### E. Rebuild and verify each layer

After the utility and parser guards, build serially to make diagnostics easier to
read. After the field-role migration and boundary-condition changes, build again:

```bash
make -j1
```

Test at least these behaviors in 1D: `Join` and `InputArrayString` produce one
value; vector and matrix defaults round-trip through `ParmParse`; vector and
matrix fields each have one component but retain distinct `_x` and `_xx` names;
and constant/expression boundary conditions work independently at `xlo` and
`xhi`. Then run the existing 2D tests to verify that the dimension guards have
not changed 2D behavior.

## Runtime and model-specific issues to address

These are not all compile errors, but they can cause invalid accesses or make
particular inputs unusable after the build succeeds.

- **Geometry validation:** [`Util.cpp`](./src/Util/Util.cpp) has validation that
  accesses the y entries of `prob_lo`/`prob_hi`. Those arrays contain only one
  element in 1D, so the checks must be dimension-guarded or rewritten as loops
  over active dimensions.
- **Initial-condition and geometry builders:** [`PointList.H`](./src/IC/PointList.H)
  explicitly rejects dimensions other than 2D. Other builders with dimension-
  specific coordinate logic include `Trig`, `Trig2`, `Sphere`, `Ellipse`,
  `Notch`, `Voronoi`, `TabulatedInterface`, and BMP/PNG-based initialization.
  Determine which are used by supported 1D cases; implement a 1D definition or
  report a clear unsupported-feature error instead of allowing invalid indexing.
- **Fracture:** [`Fracture.H`](./src/Integrator/Fracture.H) explicitly aborts for
  1D and 3D. This may remain unsupported if 1D is only intended for the flame
  application, but the limitation should be explicit and tested.
- **Other dimension assumptions:** after the first compile blockers are fixed,
  rebuild all executables and review each newly exposed error. A successful
  compile of one target is not enough if the default `make` is expected to build
  all targets.

## Changes for the AP pressure-oscillation input

The current [`input.lm.ap_pressure_oscillation_1`](./input.lm.ap_pressure_oscillation_1)
is a 2D input with a quasi-1D interpretation. It has:

- `amr.n_cell = 8 40`
- two entries each for `geometry.prob_lo`, `geometry.prob_hi`, and
  `geometry.is_periodic`
- an x-periodic, y-nonperiodic domain
- initial-condition expressions that use `y`

If the intended 1D direction is the existing long, nonperiodic y direction, map
that physical direction to 1D's only coordinate, x. That means one mesh count
(initially 40 to retain the existing axial cell count), one coordinate interval
(initially the existing y interval, `-1 mm` to `+1 mm`), and one nonperiodic
flag. Expressions that vary with the old y coordinate must use x in the 1D
input. Recheck all boundary-condition entries against the 1D BC implementation.

Do not blindly shorten every multi-value input line: many arrays describe
species, reactions, or material components rather than spatial dimensions.
Determine each field's meaning from its input schema/code. Use a new input file
and a unique `plot_file`; do not attempt to restart a 1D run from a 2D checkpoint.

This coordinate mapping is conditional on y being the intended physical
direction. Confirm the axis and boundary conditions before treating a converted
input as equivalent to the original 2D simulation.

A separate smoke-test copy,
[`input.lm.ap_pressure_oscillation_1d_smoke`](./input.lm.ap_pressure_oscillation_1d_smoke),
was created with the long, nonperiodic y direction mapped to x. It uses 40 cells,
a single nonperiodic axis, a short stop time, and a unique plot path. The original
2D input was left unchanged.

## Tests and tooling

- ALAMO's documented regression workflow is `make test`; its documented
  configure/test examples cover 2D and 3D ([README](./README.rst)). Add explicit
  1D test sections for supported executables and models.
- The test runner accepts a dimension and constructs dimension-suffixed
  executable names ([runtests.py](./scripts/runtests.py)), but the existing test
  configurations inspected here do not declare `dim=1`.
- A test named “1D” does not necessarily test a 1D executable. For example,
  [`FlowRiemannShockTube/input`](./tests/FlowRiemannShockTube/input) declares
  `dim=2` and a two-dimensional mesh; it is a quasi-1D problem, not a 1D-build
  regression test.
- Add small tests for 1D geometry parsing, x-coordinate initial conditions,
  xlo/xhi Dirichlet and Neumann conditions, periodicity, and at least one
  supported physics path. Compare a result to an analytic or trusted reference.
- The smoke copy completed two steps using
  `./bin/alamo-1d-clang++ input.lm.ap_pressure_oscillation_1d_smoke allow_unused=1 thermal.on=0 plot_file=results/phase1_pure_ap_1d_smoke/no_thermal_%DD_%Y%m%d_%H%M%S`
  and produced cell/node plotfiles under `results/phase1_pure_ap_1d_smoke/no_thermal_*`.
  This verifies startup and a short flame step only. Its `PowerLaw` model
  returns NaN for the thermal properties required by `thermal.on=1`
  ([PowerLaw.H](./src/Model/Propellant/PowerLaw.H)); many parameters from the
  legacy input are also unused by the current Flame parser. Do not interpret
  this smoke output as a validated AP result.
- To test the supported 1D thermal path, a separate
  [`input.lm.beckstead_1d_flame_surrogate`](./input.lm.beckstead_1d_flame_surrogate)
  uses `FullFeedback`, a pure-AP composition, AP thermal properties, and the
  Beckstead 500-psia/300-K ambient conditions. It ran 10,000 steps to 1 ms,
  with AMR through level 4, no NaNs, and 11 output snapshots in
  `results/beckstead_1d_flame_surrogate/output_1D_20261106_164135`.
  The eta=0.5 front moved about 0.603 um into the solid; its average speed over
  the run was about 0.603 mm/s. The front was still accelerating at the end,
  so this is not a steady burn-rate measurement.
- That case is an exploratory Flame-model surrogate, not a 1D Beckstead
  simulation: its FullFeedback kinetic/heat-fit parameters come from the
  repository's `SCPThermalSandwich` case, and its 650-K interface seed differs
  from the Beckstead initialization. The 1D target does not include the
  Beckstead LowMach/Rocfire gas chemistry. Quantitative comparison requires a
  1D LowMach/Rocfire port or a quasi-1D run with the supported 2D executable.
- Repeating the surrogate with the Beckstead 400-K ignition seed completed
  1 ms without NaNs, but did not ignite: the interface cooled below the 360-K
  cutoff, and the reported mass flux was zero by 0.5 ms. Thus the warm-seed
  run is only a numerical demonstration, not a compatible replacement for
  the Beckstead chemistry.
- The surrogate exposed a 1D AMR interpolation issue: AMReX treats nodal data
  as face-based in 1D. `NodeBilinear::interp_face` now forwards those 1D fills
  to the AMReX node interpolator. The mapper also has static lifetime because
  AMReX retains its pointer in the fill-patch coarsener cache. With both fixes,
  level-4 AMR completed the 1-ms runs, including the 400-K seed check.
- Update documentation test summaries if 1D tests are added;
  [`docs/source/Tests.py`](./docs/source/Tests.py) currently tracks 2D/3D
  indicators.
- If used to process 1D output, extend scripts that currently assume 2D/3D.
  [`scripts/analyze.py`](./scripts/analyze.py) only parses 2D/3D cell boxes, and
  [`scripts/testlib.py`](./scripts/testlib.py) creates a 2D slice for comparison.

## Recommended implementation sequence

1. Keep the existing 2D executable and input untouched. The repository supports
   dimension-specific builds; switch the active build back to 2D with:

   ```bash
   ./configure --dim=2 --comp=clang++
   make -j4
   ```

2. Resolve the 1D type-identity collisions and all scalar/elastic boundary-face
   compile errors.
3. Rebuild the 1D configuration and continue until the full intended target set
   compiles. Track each new error by root cause rather than treating repeated
   header diagnostics as separate fixes.
4. Fix runtime dimension assumptions in geometry and in the models selected for
   the intended application.
5. Add and run 1D tests before adapting the production AP input.
6. Run the converted AP case from a fresh start, inspect conservation,
   boundary behavior, flame speed, and convergence, and compare with a
   physically comparable 2D calculation.

Reconfiguring alone is not enough: the current 1D configure succeeds, but ALAMO
compilation fails before a 1D executable is linked. Avoid `make realclean` merely
to switch dimensions; it removes build artifacts and cleans the AMReX checkout.

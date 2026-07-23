# GPU-002 worked value-dispatch example

Classification: illustrative contract example. It is self-contained and not a
Hydro or Fracture implementation prescription.

## Unsafe ownership shape

```cpp
struct HostModels {
    std::vector<ModelBase*> choices;
};

// Invalid when captured state and virtual receivers were created on the host.
auto models = host_models;
amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    out(i,j,k) = models.choices[selected]->eval(in(i,j,k));
});
```

The defects are separate: a host container, host addresses, runtime virtual
dispatch, and unknown ownership/lifetime.

## Device-safe values and views

```cpp
struct TableView {
    amrex::Real const* data = nullptr;
    int size = 0;
};

struct ModelA {
    TableView coefficients;
    AMREX_GPU_HOST_DEVICE amrex::Real eval(amrex::Real x) const noexcept;
};

struct ModelB {
    amrex::Real exponent = 1.0;
    AMREX_GPU_HOST_DEVICE amrex::Real eval(amrex::Real x) const noexcept;
};

enum class ModelKind : int { a, b };

struct DeviceModels {
    amrex::GpuTuple<ModelA, ModelB> alternatives;
    ModelKind selected;

    AMREX_GPU_HOST_DEVICE
    amrex::Real eval(amrex::Real x) const noexcept {
        switch (selected) {
        case ModelKind::a: return amrex::get<0>(alternatives).eval(x);
        case ModelKind::b: return amrex::get<1>(alternatives).eval(x);
        }
        return amrex::Real(0.0); // policy reports an invalid tag separately
    }
};
```

Host parsing may still use vectors and polymorphism. Before launch it copies
variable-sized immutable coefficients into retained device storage, constructs
POD views and value alternatives, and captures only `DeviceModels`:

```cpp
amrex::Gpu::DeviceVector<amrex::Real> device_coefficients(host_coefficients.size());
amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                 host_coefficients.begin(), host_coefficients.end(),
                 device_coefficients.begin());

DeviceModels device_models{
    amrex::makeTuple(ModelA{{device_coefficients.data(), int(device_coefficients.size())}},
                     ModelB{parsed_exponent}),
    parsed_kind};

amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    out(i,j,k) = device_models.eval(in(i,j,k));
});
amrex::Gpu::streamSynchronize(); // only if required by the retained-storage boundary
```

The port must prove that all alternatives and referenced buffers remain alive,
that selection matches the CPU parser, and that no unselected alternative
silently changes behavior. A `std::variant`, switch over separate launch paths,
or another static representation is acceptable when it satisfies the same
contract. For vector-backed crack/model collections, stage a device-safe array
of values/views; do not copy the host `std::vector` object into a kernel.


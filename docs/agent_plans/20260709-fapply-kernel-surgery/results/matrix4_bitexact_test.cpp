// Standalone bit-exactness test for the Matrix4<D,Sym::Major> kernel-surgery
// edits (fapply-kernel-surgery task, Step 1):
//   * operator*(Matrix4, Set::Matrix3)  -- hand-unrolled vs original accessor loop
//   * MulCol(Matrix4, Set::Matrix, col) -- column-restricted vs operator*(...).col(col)
// Compiled once per dimension (AMREX_SPACEDIM = 2 and 3) with plain host g++.
// The REFERENCE implementations below are the ORIGINAL algorithms (accessor
// loop / full operator* + .col()); the test asserts EXACT (==) equality with
// the edited operators on randomized inputs.  Any data[] index-map or
// accumulation-order slip shows up as a bitwise mismatch.
#include <cstdio>
#include <cstdlib>
#include <random>
#include "Set/Base.H"
#include "Set/Matrix3.H"
#include "Set/Matrix4.H"
#include "Set/Matrix4_Major.H"

using Set::Sym;

// Minimal stand-ins for the ALAMO Util symbols referenced by the Matrix4
// headers, so this test links against libamrex alone (linking the real
// Util.cpp drags in the whole IO/ParmParse/Set::Constant graph, none of which
// this pure-math test exercises).  Random() drives Matrix4::Randomize();
// Abort()/globalprefix satisfy the accessor's out-of-range branch (never taken
// here).
namespace Util
{
    std::string globalprefix;
    Set::Scalar Random()
    {
        static std::mt19937_64 r(999);
        static std::uniform_real_distribution<double> d(-2.0, 2.0);
        return d(r);
    }
    void Abort(const char*) { std::abort(); }
}

// Reference: the ORIGINAL loop that operator*(Matrix4,Matrix3) used before the edit.
static Set::Vector RefMatrix3Product(const Set::Matrix4<AMREX_SPACEDIM,Sym::Major>& a,
                                     const Set::Matrix3& b)
{
    Set::Vector ret = Set::Vector::Zero();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int J = 0; J < AMREX_SPACEDIM; J++)
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                for (int L = 0; L < AMREX_SPACEDIM; L++)
                    ret(i) += a(i,J,k,L) * b(k,L,J);
    return ret;
}

int main()
{
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> dist(-3.0, 3.0);

    int fails = 0;
    const int NTRIAL = 20000;
    for (int t = 0; t < NTRIAL; ++t)
    {
        // Fill a Matrix4<Major>: use Increment() on the first trial (data[i]=i,
        // includes an exact 0.0 in data[0] -- exercises the signed-zero corner),
        // random distinct values afterwards.
        Set::Matrix4<AMREX_SPACEDIM,Sym::Major> A =
            (t == 0) ? Set::Matrix4<AMREX_SPACEDIM,Sym::Major>::Increment()
                     : Set::Matrix4<AMREX_SPACEDIM,Sym::Major>();
        if (t != 0)
        {
            // distinct random entries via the public accessor's underlying storage
            // (Randomize fills data[] directly)
            A.Randomize();
        }

        // Matrix3 with distinct values (includes negatives).
        Set::Matrix3 B3;
        for (int k = 0; k < AMREX_SPACEDIM; ++k)
            for (int l = 0; l < AMREX_SPACEDIM; ++l)
                for (int J = 0; J < AMREX_SPACEDIM; ++J)
                    B3(k,l,J) = dist(rng);

        // Matrix (for MulCol) with distinct values.
        Set::Matrix B2;
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            for (int j = 0; j < AMREX_SPACEDIM; ++j)
                B2(i,j) = dist(rng);

        // --- Test 1: Matrix4 * Matrix3 (edited operator vs original loop) ---
        Set::Vector got3  = A * B3;
        Set::Vector ref3  = RefMatrix3Product(A, B3);
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            if (got3(i) != ref3(i)) { ++fails;
                std::printf("FAIL Matrix3 trial=%d i=%d got=%.17g ref=%.17g\n", t, i, got3(i), ref3(i)); }

        // --- Test 2: MulCol vs full operator*(...).col(col) (reference) ---
        Set::Matrix full = A * B2;
        for (int c = 0; c < AMREX_SPACEDIM; ++c)
        {
            Set::Vector got = Set::MulCol(A, B2, c);
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
                if (got(i) != full(i,c)) { ++fails;
                    std::printf("FAIL MulCol trial=%d col=%d i=%d got=%.17g ref=%.17g\n", t, c, i, got(i), full(i,c)); }
        }
    }

    std::printf("[DIM=%d] %d trials, %d mismatches -> %s\n",
                AMREX_SPACEDIM, NTRIAL, fails, fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

#include <algorithm>
#include <array>
#include <cmath>
#include <stdlib.h>

#include "Set/Matrix4.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"

#include "Test/Numeric/Stencil.H"
#include "Test/Set/Matrix4.H"

#include "Operator/Elastic.H"

#include "Numeric/Interpolator/Test.H"
#include "Numeric/Interpolator/Linear.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Linear/Laplacian.H"
#include "Model/Solid/Linear/Transverse.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "Model/Solid/Finite/PseudoLinear/Cubic.H"
#include "Model/Solid/Finite/PseudoAffine/Cubic.H"
#include "Model/Solid/Linear/Hexagonal.H"
#include "Model/Solid/Affine/Hexagonal.H"
#include "Model/Chemistry/Rocfire.H"

#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/ODE/BackwardEuler.H"
#include "Solver/Local/ODE/ForwardEuler.H"
#include "Solver/Local/ODE/ODE.H"

#include "Unit/Test.H"

namespace
{
struct RocfireJacobianTestGas
{
    Set::Scalar Rg = Set::Constant::Rg;
    std::array<Set::Scalar, 6> MW = {{26.0, 28.0, 24.0, 30.0, 22.0, 32.0}};

    Set::Scalar cp_mol_species(Set::Scalar temperature, int species) const
    {
        return MW[species] *
            (800.0 + 20.0 * species + 0.02 * (species + 1) * temperature);
    }

    Set::Scalar dcp_mol_species_dT(Set::Scalar /*temperature*/, int species) const
    {
        return MW[species] * 0.02 * (species + 1);
    }
};
}

int main (int argc, char* argv[])
{
    Util::Initialize(argc, argv);

    int failed = 0;

    Util::globalprefix = "  │  ";


    #define MODELTEST(TYPE) \
        Util::Test::Message(#TYPE); \
        { \
            int subfailed = 0; \
            subfailed += Util::Test::SubMessage("PODTest",         TYPE::PODTest<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("ArithmeticTest",  TYPE::ArithmeticTest<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("DerivativeTest1", TYPE::DerivativeTest1<TYPE>(true)); \
            subfailed += Util::Test::SubMessage("DerivativeTest2", TYPE::DerivativeTest2<TYPE>(true)); \
            if (TYPE::kinvar == Model::Solid::KinematicVariable::F) \
            { \
                subfailed += Util::Test::SubMessage("MaterialFrameIndifference", TYPE::MaterialFrameIndifference<TYPE>(true)); \
            } \
            failed += Util::Test::SubFinalMessage(subfailed); \
        }
    MODELTEST(Model::Solid::Linear::Isotropic);
    MODELTEST(Model::Solid::Linear::Cubic);
    MODELTEST(Model::Solid::Linear::Laplacian);
    MODELTEST(Model::Solid::Linear::Transverse);
    MODELTEST(Model::Solid::Affine::Isotropic);
    MODELTEST(Model::Solid::Affine::Cubic);
    MODELTEST(Model::Solid::Linear::Hexagonal);
    MODELTEST(Model::Solid::Affine::Hexagonal);
    MODELTEST(Model::Solid::Finite::NeoHookean);
    MODELTEST(Model::Solid::Finite::PseudoLinear::Cubic);
    MODELTEST(Model::Solid::Finite::NeoHookeanPredeformed);
    MODELTEST(Model::Solid::Finite::PseudoAffine::Cubic);
    

    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Full>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Diagonal>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Test();
    Test::Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major>::Test();

    Util::Test::Message("Numeric::Interpolator<Linear>");
    {
        int subfailed = 0;
        Numeric::Interpolator::Test<Numeric::Interpolator::Linear<Set::Scalar> > test;
        subfailed += Util::Test::SubMessage("Match",test.Match(0));
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Numeric::Stencil test");
    {
        int subfailed = 0;
        Test::Numeric::Stencil test;
        test.Define(32);
        // first order
        subfailed += Util::Test::SubMessage("1-0-0",test.Derivative<1,0,0>(0));
        subfailed += Util::Test::SubMessage("0-1-0",test.Derivative<0,1,0>(0));
        // second order
        subfailed += Util::Test::SubMessage("2-0-0",test.Derivative<2,0,0>(0));
        subfailed += Util::Test::SubMessage("0-2-0",test.Derivative<0,2,0>(0));
        subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,2,0>(0));
        subfailed += Util::Test::SubMessage("1-1-0",test.Derivative<1,1,0>(0));
        // fourth order
        subfailed += Util::Test::SubMessage("3-1-0",test.Derivative<3,1,0>(0));
        subfailed += Util::Test::SubMessage("1-3-0",test.Derivative<1,3,0>(0));
        subfailed += Util::Test::SubMessage("2-2-0",test.Derivative<2,2,0>(0));
        subfailed += Util::Test::SubMessage("4-0-0",test.Derivative<4,0,0>(0));
        subfailed += Util::Test::SubMessage("0-4-0",test.Derivative<0,4,0>(0));
        subfailed += Util::Test::SubMessage("face gradient",
            Test::Numeric::FaceGradient());
#if AMREX_SPACEDIM>2
        // first order
        subfailed += Util::Test::SubMessage("0-0-1",test.Derivative<0,0,1>(0));
        // second order
        subfailed += Util::Test::SubMessage("0-0-2",test.Derivative<0,0,2>(0));
        subfailed += Util::Test::SubMessage("1-0-1",test.Derivative<1,0,1>(0));
        subfailed += Util::Test::SubMessage("0-1-1",test.Derivative<0,1,1>(0));
        // fourth order
        subfailed += Util::Test::SubMessage("0-0-4",test.Derivative<0,0,4>(0));
        subfailed += Util::Test::SubMessage("0-1-3",test.Derivative<0,1,3>(0));
        subfailed += Util::Test::SubMessage("0-3-1",test.Derivative<0,3,1>(0));
        subfailed += Util::Test::SubMessage("3-0-1",test.Derivative<3,0,1>(0));
        subfailed += Util::Test::SubMessage("1-0-3",test.Derivative<1,0,3>(0));
        subfailed += Util::Test::SubMessage("0-2-2",test.Derivative<0,2,2>(0));
        subfailed += Util::Test::SubMessage("2-0-2",test.Derivative<2,0,2>(0));
        subfailed += Util::Test::SubMessage("2-1-1",test.Derivative<2,1,1>(0));
        subfailed += Util::Test::SubMessage("1-2-1",test.Derivative<1,2,1>(0));
        subfailed += Util::Test::SubMessage("1-1-2",test.Derivative<1,1,2>(0));
#endif
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    //Util::Test::Message("Solver::Nonlocal::Riemann::Roe test");
    //{
    //    int subfailed = 0;
    //    subfailed += Util::Test::SubMessage("Test",Solver::Local::Riemann::Roe::Test());
    //    failed += Util::Test::SubFinalMessage(subfailed);
    //}

    Util::Test::Message("Solver::Local::ODE test");
    {
        int subfailed = 0;
        using ODESolver = Solver::Local::ODE::ODE<
            Solver::Local::ODE::ForwardEuler,
            Solver::Local::ODE::BackwardEuler>;
        ODESolver solver;
        auto rhs = [] AMREX_GPU_HOST_DEVICE(
            const Set::Scalar* state, Set::Scalar* rate)
        {
            if (!(state[0] >= 0.0)) return false;
            rate[0] = -state[0] * state[0];
            return true;
        };
        auto jacobian = [] AMREX_GPU_HOST_DEVICE(
            const Set::Scalar* state, Set::Scalar* derivative)
        {
            derivative[0] = -2.0 * state[0];
            return true;
        };

        solver.Select<Solver::Local::ODE::ForwardEuler>();
        Set::Scalar explicit_state[1] = {1.0};
        const auto explicit_result = solver.Advance<1>(
            explicit_state, 1, 0.5, rhs, jacobian, false);
        const int explicit_failed = !explicit_result.converged ||
            std::abs(explicit_state[0] - 0.5) > 1.0e-14;
        subfailed += Util::Test::SubMessage("Forward Euler", explicit_failed);

        solver.Select<Solver::Local::ODE::BackwardEuler>();
        solver.Get<Solver::Local::ODE::BackwardEuler>().Configure(
            20, 1.0e-12, 1.0e-14);
        const Set::Scalar exact = std::sqrt(3.0) - 1.0;

        Set::Scalar finite_difference_state[1] = {1.0};
        const auto finite_difference_result = solver.Advance<1>(
            finite_difference_state, 1, 0.5, rhs, jacobian, false);
        const int finite_difference_failed = !finite_difference_result.converged ||
            std::abs(finite_difference_state[0] - exact) > 1.0e-11;
        subfailed += Util::Test::SubMessage(
            "Backward Euler finite-difference Jacobian", finite_difference_failed);

        Set::Scalar analytic_state[1] = {1.0};
        const auto analytic_result = solver.Advance<1>(
            analytic_state, 1, 0.5, rhs, jacobian, true);
        const int analytic_failed = !analytic_result.converged ||
            std::abs(analytic_state[0] - exact) > 1.0e-11;
        subfailed += Util::Test::SubMessage(
            "Backward Euler analytic Jacobian", analytic_failed);
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Model::Chemistry::Rocfire Jacobian test");
    {
        int subfailed = 0;
        Model::Chemistry::Rocfire chemistry;
        chemistry.nspecies = 6;
        RocfireJacobianTestGas gas;
        constexpr int size = 6;
        const int dependent_species = Model::Chemistry::Rocfire::Primary;
        const Set::Scalar pressure = 2.0e6;
        const Set::Scalar mixture_density = 1200.0;
        const Set::Scalar gas_density = 0.6;
        const Set::Scalar temperature_scale = 1350.0;
        Set::Scalar state[size] = {0.18, 0.12, 0.15, 0.10, 0.08, 1.0};
        Set::Scalar analytic[size * size]{};

        auto rhs = [&](const Set::Scalar* z, Set::Scalar* rate)
        {
            Model::Chemistry::SpeciesArray Y{};
            Set::Scalar sum = 0.0;
            for (int column = 0; column < size - 1; ++column)
            {
                const int species = column < dependent_species ?
                    column : column + 1;
                Y[species] = z[column];
                sum += z[column];
            }
            Y[dependent_species] = 1.0 - sum;
            const Set::Scalar temperature = z[size - 1] * temperature_scale;

            Set::Scalar inverse_mw = 0.0;
            Set::Scalar cp = 0.0;
            for (int n = 0; n < size; ++n)
            {
                inverse_mw += Y[n] / gas.MW[n];
                cp += Y[n] * gas.cp_mol_species(temperature, n) / gas.MW[n];
            }
            const Set::Scalar density = pressure /
                (gas.Rg * inverse_mw * temperature);
            Model::Chemistry::SpeciesArray rhoY{};
            for (int n = 0; n < size; ++n) rhoY[n] = density * Y[n];
            const auto source = chemistry.ComputeChemistrySources(
                pressure, temperature, rhoY, 0.0, nullptr);

            for (int row = 0; row < size - 1; ++row)
            {
                const int species = row < dependent_species ? row : row + 1;
                rate[row] = source.first[species] / density;
            }
            rate[size - 1] = gas_density * source.second /
                (density * mixture_density * cp * temperature_scale);
        };

        const bool evaluated = chemistry.ComputeODEJacobian(
            state, analytic, pressure, mixture_density, gas_density, size,
            dependent_species, temperature_scale, &gas);
        Set::Scalar max_relative_error = 0.0;
        for (int column = 0; column < size && evaluated; ++column)
        {
            const Set::Scalar delta = 1.0e-6;
            Set::Scalar plus[size]{};
            Set::Scalar minus[size]{};
            Set::Scalar rate_plus[size]{};
            Set::Scalar rate_minus[size]{};
            for (int n = 0; n < size; ++n)
                plus[n] = minus[n] = state[n];
            plus[column] += delta;
            minus[column] -= delta;
            rhs(plus, rate_plus);
            rhs(minus, rate_minus);

            for (int row = 0; row < size; ++row)
            {
                const Set::Scalar numerical =
                    (rate_plus[row] - rate_minus[row]) / (2.0 * delta);
                const Set::Scalar scale = std::max(
                    1.0, std::max(std::abs(numerical),
                                  std::abs(analytic[row * size + column])));
                max_relative_error = std::max(max_relative_error,
                    std::abs(analytic[row * size + column] - numerical) / scale);
            }
        }
        const int jacobian_failed = !evaluated || max_relative_error > 1.0e-6;
        subfailed += Util::Test::SubMessage(
            "Analytic versus centered finite difference", jacobian_failed);
        failed += Util::Test::SubFinalMessage(subfailed);
    }

    Util::Test::Message("Unit test");
    {
        int subfailed = 0;
        subfailed += Util::Test::SubMessage("Equivalence", UnitTest::Equivalence(1));
        failed += subfailed;
    }


    Util::globalprefix = "";
    Util::Message(INFO,failed," tests failed");

    Util::Finalize();
    return failed;
}

#ifndef TEST_MODEL_SOLID_FINITE_NEOHOOKEANPREDEFORMED_H
#define TEST_MODEL_SOLID_FINITE_NEOHOOKEANPREDEFORMED_H

#include <array>
#include <cmath>
#include <string>

#include "Model/Solid/Finite/NeoHookeanPredeformed.H"

namespace Test
{
namespace Model
{
namespace Solid
{
namespace Finite
{
inline int
NeoHookeanPredeformedDerivatives(bool verbose = false)
{
    ::Model::Solid::Finite::NeoHookeanPredeformed model;
    model.mu = 4.0;
    model.kappa = 7.0;
    model.F0 = ::Set::Matrix::Identity();
    model.F0(0, 0) = 1.05;
    model.F0(0, 1) = 0.03;
    model.F0(1, 0) = -0.02;
    model.F0(1, 1) = 0.97;
#if AMREX_SPACEDIM == 3
    model.F0(0, 2) = -0.01;
    model.F0(1, 2) = 0.04;
    model.F0(2, 0) = 0.01;
    model.F0(2, 1) = -0.03;
    model.F0(2, 2) = 1.08;
#endif

    ::Set::Matrix F = ::Set::Matrix::Identity();
    F(0, 0) = 1.02;
    F(0, 1) = -0.04;
    F(1, 0) = 0.05;
    F(1, 1) = 0.99;
#if AMREX_SPACEDIM == 3
    F(0, 2) = 0.02;
    F(1, 2) = -0.01;
    F(2, 0) = 0.03;
    F(2, 1) = 0.01;
    F(2, 2) = 1.04;
#endif

    constexpr ::Set::Scalar dx = 1.0e-7;
    ::Set::Matrix dw_numeric = ::Set::Matrix::Zero();
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
        {
            ::Set::Matrix dF = ::Set::Matrix::Zero();
            dF(i, j) = dx;
            dw_numeric(i, j) =
                (model.W(F + dF) - model.W(F - dF)) / (2.0 * dx);
        }

    const ::Set::Matrix dw_exact = model.DW(F);
    const ::Set::Scalar dw_relative_error =
        (dw_numeric - dw_exact).norm() / dw_numeric.norm();

    ::Set::Matrix4<AMREX_SPACEDIM, ::Set::Sym::Major> ddw_numeric =
        ::Set::Matrix4<AMREX_SPACEDIM, ::Set::Sym::Major>::Zero();
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
            for (int k = 0; k < AMREX_SPACEDIM; ++k)
                for (int l = 0; l < AMREX_SPACEDIM; ++l)
                {
                    ::Set::Matrix dF = ::Set::Matrix::Zero();
                    dF(k, l) = dx;
                    ddw_numeric(i, j, k, l) =
                        (model.DW(F + dF) - model.DW(F - dF))(i, j)
                        / (2.0 * dx);
                }

    const auto ddw_exact = model.DDW(F);
    const ::Set::Scalar ddw_relative_error =
        (ddw_numeric - ddw_exact).Norm() / ddw_numeric.Norm();

    if (verbose)
    {
        Util::Message(INFO, "nonidentity F0 DW relative error = ",
                      dw_relative_error);
        Util::Message(INFO, "nonidentity F0 DDW relative error = ",
                      ddw_relative_error);
    }

    return (dw_relative_error > 1.0e-6
            || !std::isfinite(dw_relative_error))
         + (ddw_relative_error > 1.0e-4
            || !std::isfinite(ddw_relative_error));
}

inline int
NeoHookeanPredeformedFreeExpansion(bool verbose = false)
{
    ::Model::Solid::Finite::NeoHookeanPredeformed model;
    model.mu = 4.0e9;
    model.kappa = 7.0e9;
    constexpr ::Set::Scalar eigenstrain = 2.5e-3;
    model.F0 =
        (1.0 + eigenstrain) * ::Set::Matrix::Identity();

    const ::Set::Matrix stress = model.DW(model.F0);
    const ::Set::Scalar youngs_modulus =
        9.0 * model.kappa * model.mu
        / (3.0 * model.kappa + model.mu);
    const ::Set::Scalar normalized_stress =
        stress.lpNorm<Eigen::Infinity>()
        / (youngs_modulus * std::abs(eigenstrain));

    if (verbose)
        Util::Message(INFO, "free-expansion normalized stress = ",
                      normalized_stress);
    return normalized_stress <= 1.0e-6 ? 0 : 1;
}

inline int
NeoHookeanPredeformedFieldNames(bool verbose = false)
{
    ::Set::Field<::Model::Solid::Finite::NeoHookeanPredeformed> field;
    field.name = "model";
#if AMREX_SPACEDIM == 2
    const std::array<std::string, 6> expected = {
        "model_mu", "model_kappa", "model_F0xx",
        "model_F0xy", "model_F0yx", "model_F0yy"};
#elif AMREX_SPACEDIM == 3
    const std::array<std::string, 11> expected = {
        "model_mu", "model_kappa", "model_F0xx", "model_F0xy",
        "model_F0xz", "model_F0yx", "model_F0yy", "model_F0yz",
        "model_F0zx", "model_F0zy", "model_F0zz"};
#endif
    int failed = field.NComp() != static_cast<int>(expected.size());
    for (int i = 0; i < field.NComp(); ++i)
        failed += field.Name(i) != expected[i];
    if (verbose && failed)
        Util::Message(INFO, "F0 field component naming failures = ", failed);
    return failed;
}
}
}
}
}

#endif

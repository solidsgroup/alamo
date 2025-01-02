#include "IO/ParmParse.H"
#include "Integrator/ScimitarX.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "BC/Constant.H"
#include "IC/SodShock.H"
#include "Numeric/Stencil.H"
#include "Util/Util.H"
#include "Util/ScimitarX_Util.H"

namespace Integrator
{

// Define the static member variable
Util::ScimitarX_Util::getVariableIndex ScimitarX::variableIndex;

ScimitarX::ScimitarX(IO::ParmParse& pp) : ScimitarX()
{
    pp.queryclass(*this);
}

void
ScimitarX::Parse(ScimitarX& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::ScimitarX::Parse()");
    {
        ScimitarX::SolverTypeManager& setIndex = ScimitarX::SolverTypeManager::getInstance();

        std::string solverTypeStr;

        if (pp.query("SolverType", solverTypeStr)) 
        {
            auto it = ScimitarX::stringToSolverType.find(solverTypeStr);
            if (it != ScimitarX::stringToSolverType.end()) {
                value.solverType = it->second;

                ScimitarX::variableIndex = setIndex.computeAndAssignVariableIndices(value.solverType);
                int number_of_components = ScimitarX::variableIndex.NVAR_MAX;

                // Call the setupBoundaryConditions function
                IO::ParmParse bc_pp = setupPVecBoundaryConditions(pp, ScimitarX::variableIndex);

                // Example: Print consolidated boundary conditions
                std::string result;
                if (bc_pp.query("bc.pvec.type.xlo", result)) {
                    std::cout << "Boundary condition (type.xlo): " << result << std::endl;
                }
                if (bc_pp.query("bc.pvec.val.xlo", result)) {
                    std::cout << "Boundary condition (val.xlo): " << result << std::endl;
                }

                value.bc_PVec = new BC::Constant(number_of_components, pp, "bc.pvec");
                value.bc_Pressure = new BC::Constant(1, pp, "bc.pressure");
            } else {
                Util::Abort(__FILE__, __func__, __LINE__, "Invalid SolverType: " + solverTypeStr);
            }
        }
    }
    
    // Register New Fabs
    {
        value.RegisterNewFab(value.QVec_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "QVec", false);
        value.RegisterNewFab(value.QVec_old_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "QVec_old", false);

#if AMREX_SPACEDIM >= 1
        value.RegisterNewFab(value.XFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "xflux", false);
#endif
#if AMREX_SPACEDIM >= 2
        value.RegisterNewFab(value.YFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "yflux", false);
#endif
#if AMREX_SPACEDIM == 3
        value.RegisterNewFab(value.ZFlux_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "zflux", false);
#endif
        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "SourceVec", false);
        value.RegisterNewFab(value.PVec_mf, &value.bc_nothing, value.number_of_components, value.number_of_ghost_cells, "PrimitiveVec", true); 
        value.RegisterNewFab(value.Pressure_mf, &value.bc_nothing, 1, value.number_of_ghost_cells, "Pressure", true);
    }
    // Initial Conditions
    {
        std::string type = "constant";
        pp.query("ic.pvec.type", type);

        if (type == "sodshock") {
            value.ic_PVec = new IC::SodShock(value.geom, pp, "ic.pvec.sodshock", ScimitarX::variableIndex);
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid ic.pvec.type: " + type);
        }

        pp.query("ic.pressure.type", type);
        if (type == "sodshock") {
            value.ic_Pressure = new IC::SodShock(value.geom, pp, "ic.pressure.sodshock");
        } else {
            Util::Abort(__FILE__, __func__, __LINE__, "Invalid ic.pressure.type: " + type);
        }
    } 
}

void ScimitarX::Initialize(int lev)
{
    ic_PVec->Initialize(lev, PVec_mf);
    ic_Pressure->Initialize(lev, Pressure_mf);

    QVec_mf[lev]->setVal(0.0);
    QVec_old_mf[lev]->setVal(0.0);
    XFlux_mf[lev]->setVal(0.0);
#if AMREX_SPACEDIM >= 2
    YFlux_mf[lev]->setVal(0.0);
#endif
#if AMREX_SPACEDIM == 3
    ZFlux_mf[lev]->setVal(0.0);
#endif
    Source_mf[lev]->setVal(0.0);
}

void ScimitarX::Advance(int lev, Set::Scalar /*time*/, Set::Scalar /*dt*/)
{
    std::swap(*QVec_mf[lev], *QVec_old_mf[lev]);

    for (amrex::MFIter mfi(*QVec_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();

        amrex::Array4<const Set::Scalar> const& QVec_old = (*QVec_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& QVec     = (*QVec_mf[lev]).array(mfi);

        for (int n = 0; n < number_of_components; ++n) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                QVec(i, j, k, n) = QVec_old(i, j, k, n);
            });
        }
    }
}

void ScimitarX::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));
    Set::Scalar refinement_threshold = 10.0; // Set refinement threshold to 10

    for (amrex::MFIter mfi(*Pressure_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<Set::Scalar> const& pressure = (*Pressure_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=](int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(pressure, i, j, k, 0, DX);
            Set::Scalar grad_magnitude = grad.lpNorm<2>();

            if (grad_magnitude * dr > refinement_threshold) {
                tags(i, j, k) = amrex::TagBox::SET;
            }
        });
    }

    amrex::Print() << "Refinement threshold set to: " << refinement_threshold << "\n"; 
}

void ScimitarX::TimeStepBegin(Set::Scalar a_time, int a_iter) {
    // Placeholder implementation
}

void ScimitarX::TimeStepComplete(Set::Scalar time, int lev) {
    // Placeholder implementation
}

void ScimitarX::Regrid(int lev, Set::Scalar time) {
    // Placeholder implementation
}

IO::ParmParse ScimitarX::setupPVecBoundaryConditions(IO::ParmParse& pp, const Util::ScimitarX_Util::getVariableIndex& variableIndex)
{
    int n_components = variableIndex.NVAR_MAX;
    std::vector<std::string> component_names(n_components);

    for (const auto& [variable, index] : variableIndex.variableIndexMap) {
        component_names[index] = std::to_string(variable);
    }

    IO::ParmParse bc_pp;

#if AMREX_SPACEDIM == 1
    const std::vector<std::string> sides = {"xlo", "xhi"};
#elif AMREX_SPACEDIM == 2
    const std::vector<std::string> sides = {"xlo", "xhi", "ylo", "yhi"};
#elif AMREX_SPACEDIM == 3
    const std::vector<std::string> sides = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};
#else
#error "Unsupported AMREX_SPACEDIM"
#endif

    for (const std::string& side : sides) {
        std::string type_combined, val_combined;

        for (int i = 0; i < n_components; ++i) {
            std::string type_key = "bc.pvec." + component_names[i] + ".type." + side;
            std::string val_key = "bc.pvec." + component_names[i] + ".val." + side;

            std::string type_value;
            if (pp.query(type_key.c_str(), type_value)) {
                type_combined += type_value + " ";
            } else {
                type_combined += "neumann ";
            }

            std::string value_str;
            if (pp.query(val_key.c_str(), value_str)) {
                val_combined += value_str + " ";
            } else {
                val_combined += "0.0 ";
            }
        }

        type_combined = type_combined.substr(0, type_combined.size() - 1);
        val_combined = val_combined.substr(0, val_combined.size() - 1);

        bc_pp.add(("bc.pvec.type." + side).c_str(), type_combined);
        bc_pp.add(("bc.pvec.val." + side).c_str(), val_combined);
    }

    return bc_pp;
}

} // namespace Integrator


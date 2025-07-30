#ifndef IC_LS_EXPRESSION_cpp_
#define IC_LS_EXPRESSION_cpp_
#include "IC/LS/Expression.H"
#include "Integrator/NarrowBandLevelset.H" 

using namespace IC::LS;

Expression::Expression(amrex::Vector<amrex::Geometry>& _geom)
    : IC<Set::Scalar>(_geom) {}

Expression::Expression(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
    : IC<Set::Scalar>(_geom) {
    pp_queryclass(name, *this);
}

void Expression::SetTubeField(Set::Field<Set::Scalar>* tube_field_ptr) {
    Tube_mf = tube_field_ptr;
}

void Expression::SetZeroField(Set::Field<Set::Scalar>* zero_field_ptr) {
    Zero_mf = zero_field_ptr;
}

void Expression::SetCPTField(Set::Field<Set::Scalar>* cpt_field_ptr) {
    CPT_mf = cpt_field_ptr;
}

void Expression::SetLSoldField(Set::Field<Set::Scalar>* LSold_field_ptr) {
    LSold_mf = LSold_field_ptr;
}

void Expression::SetLSData(LevelSetData* data_ptr) {
    LSData = data_ptr;
}

void Expression::Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar a_time)
{
    // Define geometry constants
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
    const Set::Scalar INNERTUBE = INNERBAND_Width * min_DX;
    const Set::Scalar OUTERTUBE = BAND_WIDTH * min_DX;
    const auto& physical_domain = geom[lev].Domain();

    // Define Tube values
    const int INTERFACE   = NarrowBandTubeType::INTERFACE;
    const int INNERBAND   = NarrowBandTubeType::INNERBAND;
    const int INSIDETUBE  = NarrowBandTubeType::INSIDETUBE;
    const int EDGEPOINT   = NarrowBandTubeType::EDGEPOINT;
    const int OUTSIDETUBE = NarrowBandTubeType::OUTSIDETUBE;   

    // Define Neighbor structure
    const int num_neighbors = Neighbors::num_neighbors;
    const auto& offsets = Neighbors::offsets;

    // Define tmp field to store individual levelsets
    const auto& ba = a_field[lev]->boxArray();
    const auto& dm = a_field[lev]->DistributionMap();
    const int ncomp = f.size();
    const int nghost = a_field[lev]->nGrow();
    amrex::MultiFab tmp_field(ba, dm, ncomp, nghost);

    // Define tmp_Tube to store narrowband values for each component
    amrex::iMultiFab tmpTube_imf(ba, dm, ncomp, nghost);

    // Define vector for min/max IntVects containing narrowband cells
    amrex::Gpu::ManagedVector<amrex::IntVect> min_iv(ncomp, amrex::IntVect::TheMaxVector());
    amrex::Gpu::ManagedVector<amrex::IntVect> max_iv(ncomp, amrex::IntVect::TheMinVector());

    auto p_min_iv = min_iv.data();
    auto p_max_iv = max_iv.data();

    if (LSData) {
        if (LSData->num_objs != ncomp) LSData->num_objs = ncomp;
        LSData->objects[lev].resize(ncomp);
    }

    // Assert that number of levelset fields is equal to number max levelset ID and that length of 
    // levelset vector is equal to ncomp
    AMREX_ALWAYS_ASSERT(levelset_id.size() == ncomp);
    AMREX_ALWAYS_ASSERT(levelset_id.minCoeff() == 0);
    AMREX_ALWAYS_ASSERT(levelset_id.maxCoeff() == a_field[lev]->nComp()-1);
    
    for (amrex::MFIter mfi(tmp_field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& ghost_bx = mfi.growntilebox();

        auto tmp = tmp_field.array(mfi);
        amrex::IndexType type = tmp_field.ixType();

        auto zero_arr = (*(*Zero_mf)[lev]).array(mfi);
        auto cpt_arr = (*(*CPT_mf)[lev]).array(mfi);
        auto tmp_tube_arr = tmpTube_imf.array(mfi);

        for (int n = 0; n < ncomp; ++n) {

            // Get levelset id
            const int ils = static_cast<int>(levelset_id[n]);

            // ------------------------------------------------------------------
            // Pass 1: Compute level set values and initial tube classification
            // ------------------------------------------------------------------
            amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, IC<Set::Scalar>::geom[lev], type);
                Set::Scalar LScell = 0.0;

    #if AMREX_SPACEDIM == 1
                LScell = f[n](x(0), 0.0, 0.0, a_time);
    #elif AMREX_SPACEDIM == 2
                LScell = f[n](x(0), x(1), 0.0, a_time);
    #elif AMREX_SPACEDIM == 3
                LScell = f[n](x(0), x(1), x(2), a_time);
    #endif

                // Clamp the levelset value
                LScell = std::clamp(LScell, -OUTERTUBE, OUTERTUBE);
                tmp(i,j,k,n) = LScell;

                // Set the CPT flag
                if (LScell <= 0.0) cpt_arr(i,j,k,ils) = n + 1;

                // Get object narrowband value and initialize to INSIDE/OUTSIDE
                tmp_tube_arr(i,j,k,n) = (std::abs(LScell) < OUTERTUBE) ? INSIDETUBE : OUTSIDETUBE;  
            });

            // ------------------------------------------------------------------
            // Pass 2: Tag tube types for object
            // ------------------------------------------------------------------
            amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar LScell   = tmp(i,j,k,n);
                Set::Scalar absG     = std::abs(LScell);
                int& tube_val = tmp_tube_arr(i,j,k,n);

                // Skip OUTSIDE cells
                if (tube_val == OUTSIDETUBE) return;

                // Step 1: Tag INTERFACE or INNERBAND
                if (absG < INNERTUBE) {
                    bool is_interface = false;

                    for (int d = 0; d < num_neighbors; ++d) {
                        int ni = i + offsets[d][0];
                        int nj = j + offsets[d][1];
                        int nk = k + offsets[d][2];

                        if (!ghost_bx.contains(ni,nj,nk)) continue;

                        Set::Scalar LSnbr = tmp(ni,nj,nk,n);
                        if (LScell * LSnbr <= 0.0) {
                            tmp_tube_arr(i,j,k,n) = INTERFACE;
                            zero_arr(i,j,k,ils) = n + 1;
                            is_interface = true;
                            break;
                        }
                    }

                    if (!is_interface) tube_val = INNERBAND;
                }

                // Step 2: Tag EDGEPOINTs (self)
                if (tube_val == INSIDETUBE) {
                    for (int d = 0; d < num_neighbors; ++d) {
                        int ni = i + offsets[d][0];
                        int nj = j + offsets[d][1];
                        int nk = k + offsets[d][2];

                        if (!ghost_bx.contains(ni,nj,nk)) continue;

                        if (tmp_tube_arr(ni,nj,nk,n) == OUTSIDETUBE) {
                            tmp_tube_arr(ni,nj,nk,n) = EDGEPOINT;
                        }
                    }
                }
            });  
            
            // ----------------------------------
            // Pass 3: Bounding box for EDGEPOINTs
            // ----------------------------------
            amrex::ParallelFor(ghost_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (tmp_tube_arr(i,j,k,n) == EDGEPOINT) {
                    amrex::IntVect coord(AMREX_D_DECL(i,j,k));
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        amrex::Gpu::Atomic::Min(&(p_min_iv[n][d]), coord[d]);
                        amrex::Gpu::Atomic::Max(&(p_max_iv[n][d]), coord[d]);
                    }
                }
            });
        }
    }

    // Fill boundaries
    tmp_field.FillBoundary();
    tmpTube_imf.FillBoundary();
    (*(*Zero_mf)[lev]).FillBoundary();
    (*(*CPT_mf)[lev]).FillBoundary();

    // 1. Group component indices by their levelset_id
    std::map<int, std::vector<int>> ils_to_comps;
    for (int n = 0; n < ncomp; ++n) {
        ils_to_comps[static_cast<int>(levelset_id[n])].push_back(n);
    }

    // 2. For each group, compute min over components and write to LS(i,j,k,ils)
    for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        auto tmp = tmp_field.const_array(mfi);
        auto LS = a_field.Patch(lev, mfi);
        auto LSold = LSold_mf->Patch(lev, mfi);

        // Iterate over groups
        for (const auto& [ils, comp_list] : ils_to_comps) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Scalar val = tmp(i,j,k,comp_list[0]);  // first comp in group
                for (int idx = 1; idx < comp_list.size(); ++idx) {
                    val = std::min(val, tmp(i,j,k,comp_list[idx]));
                }
                LS(i,j,k,ils) = val;
                LSold(i,j,k,ils) = val;
            });
        }
    }

    a_field[lev]->FillBoundary(); 
    (*(*LSold_mf)[lev]).FillBoundary();

    // Min over Tube
    for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();

        const auto& tmp_arr = tmpTube_imf.const_array(mfi);
        const auto& tube_arr = (*(*Tube_mf)[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int n = 0; n < ncomp; ++n) {
                // Get levelset id
                const int ils = static_cast<int>(levelset_id[n]);

                int Tube_val = tmp_arr(i,j,k,n);
                amrex::Gpu::Atomic::Min(&(tube_arr(i,j,k,ils)), static_cast<Set::Scalar>(Tube_val));
            }
        });
    }

    (*(*Tube_mf)[lev]).FillBoundary();

    // Update LSData.objs
    for (int n = 0; n < ncomp; ++n) {
        auto& obj = LSData->objects[lev][n];
        // Define obj meta
        obj.FlowData = LSData;
        obj.object_id = n + 1;
        obj.ls_id = static_cast<int>(levelset_id[n]);
        obj.has_narrowband = (min_iv[n] <= max_iv[n]);

        // Define obj geometry
        amrex::Box Tube_domain(min_iv[n], max_iv[n]);
        // Clip to physical domain
        Tube_domain = Tube_domain & physical_domain;
        obj.Tube_domain = Tube_domain;

        // Allocate and initialize fields
        obj.AllocateFabs(tmpTube_imf, n);
        obj.LoadLSFields(*a_field[lev], *a_field[lev]);
        obj.LoadIntFields((*(*Tube_mf)[lev]), (*(*Zero_mf)[lev]), (*(*CPT_mf)[lev]));
    }
}

void IC::LS::Expression::Parse(Expression& value, IO::ParmParse& pp) {
    std::string coordstr = "";
    // coordinate system to use
    pp_query_validate("coord", coordstr, {"cartesian","polar"}); 
    if (coordstr == "cartesian") value.coord = Expression::CoordSys::Cartesian;
    else if (coordstr == "polar") value.coord = Expression::CoordSys::Polar;
    else Util::Abort(INFO, "unsupported coordinates ", coordstr);

    // Number of levelsets
    pp.queryarr("levelset_id", value.levelset_id);

    std::vector<std::string> expression_strs;
    // Mathematical expression in terms of x,y,z,t (if coord=cartesian)
    // or r,theta,z,t (if coord=polar) and any defined constants.
    pp.query_enumerate("region", expression_strs);

    for (unsigned int i = 0; i < expression_strs.size(); i++)
    {
        value.parser.push_back(amrex::Parser(expression_strs[i]));

        //
        // Read in user-defined constants and add them to the parser
        //
        std::string prefix = pp.getPrefix();
        std::set<std::string> entries = pp.getEntries(prefix + ".constant");//"constant");
        std::set<std::string>::iterator entry;
        for (entry = entries.begin(); entry != entries.end(); entry++)
        {
            IO::ParmParse pp;
            std::string fullname = *entry;
            Set::Scalar val  = NAN;
            pp_query(fullname.data(),val);
            std::string name = Util::String::Split(fullname,'.').back();
            value.parser.back().setConstant(name,val);
        }

        if (value.coord == Expression::CoordSys::Cartesian)
        {
            value.parser.back().registerVariables({ "x","y","z","t" });
            value.f.push_back(value.parser.back().compile<4>());
        }
        else if (value.coord == Expression::CoordSys::Polar)
        {
            value.parser.back().registerVariables({ "r","theta","z","t" });
            value.f.push_back(value.parser.back().compile<4>());
        }
    }
}

#endif

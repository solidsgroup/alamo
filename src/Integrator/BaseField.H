#ifndef INTEGRATOR_BASEFIELD_H
#define INTEGRATOR_BASEFIELD_H

#include "Numeric/Interpolator/NodeBilinear.H"

namespace Integrator
{

class BaseField
{
public:
    virtual void RemakeLevel (int lev,       
                            amrex::Real time, 
                            const amrex::BoxArray& cgrids, 
                            const amrex::DistributionMapping& dm) = 0;
    virtual void MakeNewLevelFromCoarse (int lev, 
                                        amrex::Real time, 
                                        const amrex::BoxArray& cgrids, 
                                        const amrex::DistributionMapping& dm) = 0;
    virtual void MakeNewLevelFromScratch (int lev, 
                                        amrex::Real t, 
                                        const amrex::BoxArray& cgrids,
                                        const amrex::DistributionMapping& dm) = 0;
    virtual void SetFinestLevel(const int a_finestlevel) = 0;
                      
};

template<class T>
class Field : public BaseField
{
public:
    class EmptyBC
    {
    public:
        virtual void operator () (amrex::FabArray<amrex::BaseFab<T>> & /*mf*/, 
                                int /*dcomp*/, int /*ncomp*/, amrex::IntVect const& /*nghost*/, 
                                amrex::Real /*time*/, int /*bccomp*/) 
        {  /* Do nothing - this is a shell to satisfy the FillPatch functions.*/   }
        amrex::BCRec GetBCRec() {return amrex::BCRec();}
    };

    Field(Set::Field<T> & a_field, 
        const amrex::Vector<amrex::Geometry> &a_geom, 
        const amrex::Vector<amrex::IntVect> &a_refRatio,
        int a_ncomp, int a_nghost) : 
        m_field(a_field), m_geom(a_geom), m_refRatio(a_refRatio), 
        m_ncomp(a_ncomp), m_nghost(a_nghost)
    {} 
    
    void 
    FillPatch (int lev, amrex::Real time,
                amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>> &source_mf,
                amrex::FabArray<amrex::BaseFab<T>> &destination_mf, 
                int icomp)
    {
        EmptyBC physbc;
        if (lev == 0)
        {
            amrex::Vector<amrex::FabArray<amrex::BaseFab<T>>*> smf;
            smf.push_back(source_mf[lev].get());
            amrex::Vector<amrex::Real> stime;
            stime.push_back(time);
            //physbc.define(m_geom[lev]);
            amrex::FillPatchSingleLevel(destination_mf, time, smf, stime,
                                        0, icomp, destination_mf.nComp(), m_geom[lev],
                                        physbc, 0);
        } 
        else
        {
            amrex::Vector<amrex::FabArray<amrex::BaseFab<T>>*> cmf, fmf;
            cmf.push_back(source_mf[lev-1].get());
            fmf.push_back(source_mf[lev].get());
            amrex::Vector<amrex::Real> ctime, ftime;
            ctime.push_back(time);
            ftime.push_back(time);

            //physbc.define(m_geom[lev]);

            Numeric::Interpolator::NodeBilinear<T> *mapper; 

            if (destination_mf.boxArray().ixType() == amrex::IndexType::TheNodeType())
                mapper = new Numeric::Interpolator::NodeBilinear<T>();
            else
                Util::Abort(INFO,"Currently works for node-based templated fabs only!");
            //mapper = &amrex::cell_cons_interp;

            amrex::Vector<amrex::BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); 
            amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
                                    0, icomp, destination_mf.nComp(), m_geom[lev-1], m_geom[lev],
                                    physbc, 0,
                                    physbc, 0,
                                    m_refRatio[lev-1],
                                    mapper, bcs, 0);
        }

    }
    void
    FillCoarsePatch (int lev,
                    amrex::Real time,
                    int icomp,
                    int ncomp)
    {
        EmptyBC physbc;

        BL_PROFILE("Integrator::FillCoarsePatch");
        AMREX_ASSERT(lev > 0);
        amrex::Vector<amrex::FabArray<amrex::BaseFab<T>> *> cmf;
        cmf.push_back(m_field[lev-1].get());
        amrex::Vector<amrex::Real> ctime;
        ctime.push_back(time);
    
        //physbc.define(geom[lev]);
        //amrex::Interpolater* mapper = &amrex::cell_cons_interp;
        Numeric::Interpolator::NodeBilinear<T> *mapper = new Numeric::Interpolator::NodeBilinear<T>();

        amrex::Vector<amrex::BCRec> bcs(ncomp, physbc.GetBCRec());
        amrex::InterpFromCoarseLevel(*m_field[lev], time, *cmf[0], 0, icomp, ncomp, 
                                    m_geom[lev-1], m_geom[lev],
                                    physbc, 0,
                                    physbc, 0,
                                    m_refRatio[lev-1],
                                    mapper, bcs, 0);
    }


    virtual void RemakeLevel (int lev,       
                            amrex::Real time, 
                            const amrex::BoxArray& cgrids, 
                            const amrex::DistributionMapping& dm) override
    {
        if (m_field[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType())
        {
            amrex::BoxArray ngrids = cgrids;
            ngrids.convert(amrex::IntVect::TheNodeVector());
            amrex::FabArray<amrex::BaseFab<T>> new_state(ngrids, dm, m_ncomp, m_nghost);
            this->FillPatch(lev, time, m_field, new_state, 0);
            std::swap(new_state, *m_field[lev]);
        }
        else
        {
            amrex::FabArray<amrex::BaseFab<T>> new_state(cgrids, dm, m_ncomp, m_nghost);
            this->FillPatch(lev, time, m_field, new_state, 0);
            std::swap(new_state, *m_field[lev]);
        }
    }

    virtual void MakeNewLevelFromCoarse (int lev, 
                                        amrex::Real time, 
                                        const amrex::BoxArray& cgrids, 
                                        const amrex::DistributionMapping& dm) override
    {
        //if (m_field[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType())
        //{
        amrex::BoxArray ngrids = cgrids;
        ngrids.convert(amrex::IntVect::TheNodeVector());
        m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(ngrids,dm,m_ncomp,m_nghost));
        FillCoarsePatch(lev,time,0,m_ncomp);
        //}
        //else
        //{
        //    m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(cgrids,dm,m_ncomp,m_nghost));
        //    FillCoarsePatch(lev,time,0,m_ncomp);
        //}
    }

    virtual void MakeNewLevelFromScratch (int lev, 
                                        amrex::Real /*time*/, 
                                        const amrex::BoxArray& cgrids,
                                        const amrex::DistributionMapping& dm)
    {
        //if (m_field[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType())
        {
            amrex::BoxArray ngrids = cgrids;
            ngrids.convert(amrex::IntVect::TheNodeVector());
            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(ngrids,dm,m_ncomp,m_nghost));
        }
//      else
//      {
//            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(cgrids,dm,m_ncomp,m_nghost));
//      }
    }                                 

    virtual void SetFinestLevel(const int a_finestlevel)
    {
        m_field.finest_level = a_finestlevel;
    }
    
private:

    Set::Field<T> &m_field;
    const amrex::Vector<amrex::Geometry>  &m_geom;
    const amrex::Vector<amrex::IntVect> &m_refRatio;
    const int m_ncomp, m_nghost;
};

}

#endif

#ifndef INTEGRATOR_BASEFIELD_H
#define INTEGRATOR_BASEFIELD_H

#include "AMReX_FillPatchUtil.H"

#include "Util/Util.H"
#include "Set/Set.H"
#include "Numeric/Interpolator/NodeBilinear.H"

namespace Integrator
{

class BaseField
{
public:
    virtual ~BaseField() = default;

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
    virtual void FillPatch(const int lev, const Set::Scalar time)=0;
    virtual void FillBoundary(const int lev, Set::Scalar time)=0;
    virtual void AverageDown(const int lev, amrex::IntVect refRatio) = 0;

    virtual int NComp() = 0;
    virtual void Copy(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) = 0;
    bool writeout = false;
    virtual std::string Name(int) = 0;
    virtual void setName(std::string a_name) = 0;
    bool evolving = true;
    virtual void setBC(void * a_bc) = 0;
    virtual void * getBC() = 0;
    Set::Hypercube m_gridtype;
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

            if (m_gridtype == Set::Hypercube::Cell)
            {
                m_bc->define(m_geom[lev]);
                amrex::FillPatchSingleLevel(destination_mf, time, smf, stime,
                                            0, icomp, destination_mf.nComp(), m_geom[lev],
                                            *m_bc, 0);
            }
            else if (m_gridtype == Set::Hypercube::Node)
            {
                amrex::FillPatchSingleLevel(destination_mf, time, smf, stime,
                                0, icomp, destination_mf.nComp(), m_geom[lev],
                                physbc, 0);
            }
        } 
        else
        {
            amrex::Vector<amrex::FabArray<amrex::BaseFab<T>>*> cmf, fmf;
            cmf.push_back(source_mf[lev-1].get());
            fmf.push_back(source_mf[lev].get());
            amrex::Vector<amrex::Real> ctime, ftime;
            ctime.push_back(time);
            ftime.push_back(time);

            amrex::Vector<amrex::BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); 

            if (destination_mf.boxArray().ixType() == amrex::IndexType::TheNodeType())
            {
                Numeric::Interpolator::NodeBilinear<T> *mapper
                    = new Numeric::Interpolator::NodeBilinear<T>();
            
                amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
                                        0, icomp, destination_mf.nComp(), m_geom[lev-1], m_geom[lev],
                                        physbc, 0,
                                        physbc, 0,
                                        m_refRatio[lev-1],
                                        mapper, bcs, 0);

            }
            else
            {
                Util::Abort(INFO,"This implementation is not yet complete.");
                /*
                Numeric::Interpolator::CellBilinear<T> *mapper 
                    = new Numeric::Interpolator::CellBilinear<T>();
                m_bc->define(m_geom[lev]);
                amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
                                        0, icomp, destination_mf.nComp(), m_geom[lev-1], m_geom[lev],
                                        *m_bc, 0,
                                        *m_bc, 0,
                                        m_refRatio[lev-1],
                                        mapper, bcs, 0);
                */
            }
        }

    }
    void 
    FillPatch (int lev, amrex::Real time) override
    {
        this->FillPatch(lev,time,m_field,*(m_field[lev]),0);
    }
    void
    FillCoarsePatch (int lev,
                    amrex::Real time,
                    int icomp,
                    int ncomp)
    {
        BL_PROFILE("Integrator::FillCoarsePatch");
        AMREX_ASSERT(lev > 0);

        if (m_gridtype == Set::Hypercube::Node)
        {
            Util::Assert(INFO, TEST(m_field[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType()));
            EmptyBC physbc;
            amrex::Vector<amrex::FabArray<amrex::BaseFab<T>> *> cmf;
            cmf.push_back(m_field[lev-1].get());
            amrex::Vector<amrex::Real> ctime;
            ctime.push_back(time);
            Numeric::Interpolator::NodeBilinear<T> *mapper = new Numeric::Interpolator::NodeBilinear<T>();
            amrex::Vector<amrex::BCRec> bcs(ncomp, physbc.GetBCRec());
            amrex::InterpFromCoarseLevel(*m_field[lev], time, *cmf[0], 0, icomp, ncomp, 
                                        m_geom[lev-1], m_geom[lev],
                                        physbc, 0,
                                        physbc, 0,
                                        m_refRatio[lev-1],
                                        mapper, bcs, 0);
        }
        else if (m_gridtype == Set::Hypercube::Cell)
        {
            Util::Abort(INFO,"This implementation is not yet complete");
            /*
            amrex::Vector<amrex::FabArray<amrex::BaseFab<T>> *> cmf;
            cmf.push_back(m_field[lev-1].get());
            amrex::Vector<amrex::Real> ctime;
            ctime.push_back(time);
            m_bc->define(m_geom[lev]);
            Numeric::Interpolator::CellBilinear<T> *mapper = new Numeric::Interpolator::CellBilinear<T>();
            amrex::Vector<amrex::BCRec> bcs(ncomp, m_bc->GetBCRec());
            amrex::InterpFromCoarseLevel(*m_field[lev], time, *cmf[0], 0, icomp, ncomp, 
                                        m_geom[lev-1], m_geom[lev],
                                        *m_bc, 0,
                                        *m_bc, 0,
                                        m_refRatio[lev-1],
                                        mapper, bcs, 0);
            */
        }
        else
            Util::Abort(INFO,"Invalid grid");
    }


    virtual void RemakeLevel (int lev,       
                            amrex::Real time, 
                            const amrex::BoxArray& cgrids, 
                            const amrex::DistributionMapping& dm) override
    {
        if (m_gridtype == Set::Hypercube::Node)
        {
            Util::Assert(INFO, TEST(m_field[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType()));
            amrex::BoxArray ngrids = cgrids;
            ngrids.convert(amrex::IntVect::TheNodeVector());
            amrex::FabArray<amrex::BaseFab<T>> new_state(ngrids, dm, m_ncomp, m_nghost);
            this->FillPatch(lev, time, m_field, new_state, 0);
            std::swap(new_state, *m_field[lev]);
        }
        else if (m_gridtype == Set::Hypercube::Cell)
        {
            amrex::FabArray<amrex::BaseFab<T>> new_state(cgrids, dm, m_ncomp, m_nghost);
            this->FillPatch(lev, time, m_field, new_state, 0);
            std::swap(new_state, *m_field[lev]);
        }
        else
            Util::Abort(INFO,"Invalid grid");
    }

    virtual void MakeNewLevelFromCoarse (int lev, 
                                        amrex::Real time, 
                                        const amrex::BoxArray& cgrids, 
                                        const amrex::DistributionMapping& dm) override
    {
        if (m_gridtype == Set::Hypercube::Node)
        {
            amrex::BoxArray ngrids = cgrids;
            ngrids.convert(amrex::IntVect::TheNodeVector());
            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(ngrids,dm,m_ncomp,m_nghost));
            FillCoarsePatch(lev,time,0,m_ncomp);
        }
        else if (m_gridtype == Set::Hypercube::Cell)
        {
            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(cgrids,dm,m_ncomp,m_nghost));
            FillCoarsePatch(lev,time,0,m_ncomp);
        }
        else
            Util::Abort(INFO,"Invalid grid");
    }

    virtual void MakeNewLevelFromScratch (int lev, 
                                        amrex::Real /*time*/, 
                                        const amrex::BoxArray& cgrids,
                                        const amrex::DistributionMapping& dm) override
    {
        if (m_gridtype == Set::Hypercube::Node)
        {
            amrex::BoxArray ngrids = cgrids;
            ngrids.convert(amrex::IntVect::TheNodeVector());
            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(ngrids,dm,m_ncomp,m_nghost));
            m_field[lev]->setVal(T::Zero());
        }
        else if (m_gridtype == Set::Hypercube::Cell)
        {
            m_field[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(cgrids,dm,m_ncomp,m_nghost));
            m_field[lev]->setVal(T::Zero());
        }
        else
            Util::Abort(INFO,"Invalid grid");
    }                                 

    virtual void SetFinestLevel(const int a_finestlevel) override
    {
        m_field.finest_level = a_finestlevel;
    }

    virtual void FillBoundary(const int lev, Set::Scalar time) override
    {
        Util::RealFillBoundary(*m_field[lev],m_geom[lev]);
        if (m_gridtype == Set::Hypercube::Cell)
        {
            m_bc->define(m_geom[lev]);
            for (amrex::MFIter mfi(*m_field[lev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.tilebox();
                amrex::BaseFab<T>  &patch = (*(m_field[lev]))[mfi];
                m_bc->FillBoundary(patch,box,m_nghost,0,m_ncomp,time);
            }
            //
        }
    }
    
    virtual void AverageDown(const int lev, amrex::IntVect refRatio) override
    {
        if (m_gridtype == Set::Hypercube::Node)
            amrex::average_down_nodal(*m_field[lev+1],*m_field[lev],refRatio);
        else if (m_gridtype == Set::Hypercube::Cell)
            Util::Abort(INFO,"This is not yet working");
            //amrex::average_down(*m_field[lev+1],*m_field[lev],0,m_ncomp,refRatio);
            //amrex::average_down(*m_field[lev+1],*m_field[lev],0,m_ncomp,refRatio);
        else
            Util::Abort(INFO,"Incorrect grid");
    }
    
private:

    Set::Field<T> &m_field;
    const amrex::Vector<amrex::Geometry>  &m_geom;
    const amrex::Vector<amrex::IntVect> &m_refRatio;
    const int m_ncomp, m_nghost;
    BC::BC<T> *m_bc;
    
public:
    virtual int NComp() override {
        return m_field.NComp();
    }    
    virtual void Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) override
    {
        m_field.Copy(a_lev,a_dst,a_dstcomp,a_nghost);
    }
    virtual void setName(std::string a_name) override {
        m_field.name = a_name;
    }
    virtual std::string Name(int i) override {
        return m_field.Name(i);
    }
    virtual void setBC(void * a_bc) override {
        m_bc = (BC::BC<T> *)(a_bc);
    }
    virtual void * getBC() override {
        return (void*)m_bc;
    }

};

}

#endif

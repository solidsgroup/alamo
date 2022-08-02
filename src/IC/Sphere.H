#ifndef IC_SPHERE_H_
#define IC_SPHERE_H_

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"

/// \class Sphere
/// \brief Initialize a spherical inclusion
namespace IC
{
    class Sphere : public IC
    {
    public:
        /// Sphere is also used for a cylindrical case, but switching from the default type to the desired orientation.
        /// For instance, XY refers to a cylinder with the circular faces in the XY plane.
        enum Type
        {
            XYZ,
            YZ,
            ZX,
            XY
        };

        Sphere(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
        Sphere(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp) : IC(_geom)
        {pp.queryclass(*this);}
        Sphere(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : IC(_geom)
        {pp.queryclass(name,*this);}
        /// Constructor defining radius, center, dimension/orientation, and field values within and outside of the area.
        Sphere(amrex::Vector<amrex::Geometry> &_geom, Set::Scalar _radius, Set::Vector _center, Type _type = Type::XYZ, Set::Scalar _alpha_in = 1, Set::Scalar _alpha_out = 0)
            : IC(_geom)
        {
            Define(_radius, _center, _type, _alpha_in, _alpha_out);
        }

        void Define(Set::Scalar a_radius,
                    Set::Vector a_center,
                    Type a_type,
                    Set::Scalar a_alpha_in,
                    Set::Scalar a_alpha_out)
        {
            radius = a_radius;
            center = a_center;
            type = a_type;
            alpha_in = a_alpha_in;
            alpha_out = a_alpha_out;
        }

        void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
        {
            bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));
            int ncomp = a_field[0]->nComp();

            for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                // bx.grow(a_field[lev]->nGrow());

                amrex::Array4<Set::Scalar> const &field = a_field[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {

                    Set::Scalar AMREX_D_DECL(x,y,z);
                    if (cellcentered)
                    {
                        AMREX_D_TERM(x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                    y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                    z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                    }
                    else
                    {
                        AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
                                    y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
                                    z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];);
                    }

                    Set::Scalar rsq = NAN;

                    if (type == Type::XYZ)
                    {
                        // 3D Sphere
                        rsq =
                        AMREX_D_TERM(  (x-center(0))*(x-center(0)) , 
                                    + (y-center(1))*(y-center(1)) ,
                                    + (z-center(2))*(z-center(2)));
                    }
                    else if (type == Type::XY)
                    {
                        // Cylinder along Z axis
                        rsq =
                        AMREX_D_TERM(  (x-center(0))*(x-center(0)) , 
                                    + (y-center(1))*(y-center(1)) ,
                                    );
                    }
                    else if (type == Type::YZ)
                    {
                        // Cylinder along X axis
                        rsq = 
                        AMREX_D_TERM( , 
                                    + (y-center(1))*(y-center(1)) ,
                                    + (z-center(2))*(z-center(2)));
                    }
                    else if (type == Type::ZX)
                    {
                        // Cylinder along Y axis
                        rsq =
                        AMREX_D_TERM(  (x-center(0))*(x-center(0)) , 
                                    ,
                                    + (z-center(2))*(z-center(2)));
                    }

                    if ( rsq < radius*radius)
                        {
                            field(i,j,k,0) = alpha_in;     
                            if (ncomp>1) field(i,j,k,1) = 0.;     
                        }
                    else
                        {
                            field(i,j,k,0) = alpha_out;     
                            if (ncomp>1) field(i,j,k,1) = 1;     
                        } });
            }
        };
        // This is a somewhat antiquated IC that will eventually be replaced
        // with the Expression IC.
        static void Parse(Sphere &value, IO::ParmParse &pp)
        {
            pp.query("radius", value.radius); // Radius of the sphere
            pp.queryarr("center", value.center); // Vector location of the sphere center
            pp.query("inside", value.alpha_in); // Value of the field inside the sphere
            pp.query("outside", value.alpha_out); // Value of the field outside teh sphere
            std::string type;
            pp.query("type", type); // Type - can be cylinder oriented along the x, y, z directions or full sphere.
            if (type == "yz")
                value.type = Type::YZ;
            if (type == "zx")
                value.type = Type::ZX;
            if (type == "xy")
                value.type = Type::XY;
            if (type == "xyz")
                value.type = Type::XYZ;
        }

    private:
        int dim = AMREX_SPACEDIM;
        Set::Vector center = Set::Vector::Zero();
        Set::Scalar radius = 1.0;
        Set::Scalar alpha_in = 1.0;
        Set::Scalar alpha_out = 0.0;
        Type type = Type::XYZ;
    };
}
#endif

#ifndef IC_DYNAMICICNUCLEATION_H_
#define IC_DYNAMICICNUCLEATION_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "IO/ParmParse.H"

namespace IC
{
class DynamicNucleation : public IC
{
public:
    DynamicNucleation(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    // Instead of parsing from a file, we dynamically add positions and radii
    void AddNucleationPoint(const Set::Vector& position, Set::Scalar radius) {
        nucleation_positions.push_back(position); // Store the nucleation point
        radii.push_back(radius);                  // Store the corresponding radius
    }

    // Method that applies the nucleation points to the domain (i.e., adds voids/circles)
    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar eps = 1E-16) {
        Set::Vector size = Set::Size(geom[lev]);
        int ncomp = a_phi[lev]->nComp();

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx;
            amrex::IndexType type = a_phi[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unknown index type");

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                Set::Scalar min_grain_id = 0;

                for (unsigned int n = 0; n < nucleation_positions.size(); n++) {
                    Set::Scalar d = (x - nucleation_positions[n]).lpNorm<2>();  // Euclidean distance

                    Set::Scalar d1 = std::numeric_limits<Set::Scalar>::infinity();
                    if (geom[0].isPeriodic(0)) {
                        d1 = std::min((x - nucleation_positions[n] + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                                      (x - nucleation_positions[n] - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                    }
#if AMREX_SPACEDIM > 1
                    if (geom[0].isPeriodic(1)) {
                        d1 = std::min((x - nucleation_positions[n] + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                                      (x - nucleation_positions[n] - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                    }
#endif
                    d = std::min(d, d1);

                    // Add nucleation (void) based on the distance to the nucleation point
                    if (d <= (radii[n] + eps)) {
                        Set::Scalar m = 0.5 * (1 + erf((-d + radii[n]) / eps));
                        min_grain_id = min_grain_id + m * (1.0 - min_grain_id);  // Smooth addition of void
                    }
                }

                // Apply to the field (damage, phi, or other relevant field)
                phi(i, j, k, 0) = min_grain_id;
                if (ncomp > 1) phi(i, j, k, 1) = 1.0 - min_grain_id;  // Complementary field, if needed
            });
        }
    }

    void ClearNucleationPoints() {
        nucleation_positions.clear();  // Clear all existing nucleation points
        radii.clear();
    }

private:
    std::vector<Set::Vector> nucleation_positions;  // Nucleation centers (positions)
    std::vector<Set::Scalar> radii;                 // Radii of the nucleated voids
};
}
#endif  // IC_DYNAMICICNUCLEATION_H_


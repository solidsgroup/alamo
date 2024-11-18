#include "Model/Solid/Finite/Adhesion.H" // Include Adhesion header file
#include <iostream>                     // For printing outputs
#include <Eigen/Dense>                  // Include Eigen if Set::Matrix uses Eigen

int main() {
    // Create an instance of the Adhesion material
    Model::Solid::Finite::Adhesion adhesionMaterial;

    // Define material parameters
    adhesionMaterial.d = 1.0;      // Scaling factor
    adhesionMaterial.mu = 0.5;    // Shear modulus
    adhesionMaterial.kappa = 1.0; // Bulk modulus
    adhesionMaterial.zeta = 0.1;  // Penetration resistance
    adhesionMaterial.n = 2.0;     // Exponent for penetration term

    // Create a deformation gradient tensor
    Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // Assuming Set::Matrix uses Eigen
    F(0, 0) = 1.2; // Example stretch in x-direction
    F(1, 1) = 0.8; // Example compression in y-direction
    F(2, 2) = 1.0; // No deformation in z-direction

    // Compute the energy W(F)
    Set::Scalar energy = adhesionMaterial.W(F);
    std::cout << "Energy W(F): " << energy << std::endl;

    // Compute the first derivative DW(F)
    Eigen::Matrix3d derivative = adhesionMaterial.DW(F); // Assuming Set::Matrix uses Eigen
    std::cout << "First derivative DW(F):" << std::endl;
    std::cout << derivative << std::endl;

    // Compute the second derivative DDW(F) (optional)
    // If DDW(F) is a placeholder, this line might not be fully functional
    auto secondDerivative = adhesionMaterial.DDW(F);
    std::cout << "Second derivative DDW(F): (Output skipped due to placeholder)" << std::endl;

    return 0;
}

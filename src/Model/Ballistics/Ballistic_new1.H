#ifndef MODEL_BALLISTICS_BALLISTIC_NEW_H_
#define MODEL_BALLISTICS_BALLISTIC_NEW_H_

#include "IO/ParmParse.H"
//#include "Integrator/Flame.cpp"
#include <vector>      // Include vector for dynamic array
#include <cmath>       // Include cmath for sqrt

namespace Model
{
namespace Ballistics
{
    class Ballistic_new {
        public:
            Ballistic_new() {}
            virtual ~Ballistic_new() {}
        
            Set::Scalar Compute(Set::Scalar mdot, Set::Scalar pressure, Set::Scalar dt,
                                Set::Scalar At, Set::Scalar gamma, Set::Scalar R,
                                Set::Scalar T0, Set::Scalar volume) const
            {
                Set::Scalar mdot_in = mdot;
                Set::Scalar mdot_out = (pressure * At) / sqrt(gamma * R * T0);
        
                std::cout << "pressure: " << pressure << "\n";
                std::cout << "mdot_in: " << mdot_in << "\n";
                std::cout << "mdot_out: " << mdot_out << "\n";
        
                Set::Scalar dpdt = (mdot_in - mdot_out) / (volume / (R * T0));
                return pressure + dpdt * dt;
            }
        };
        
        // Chamber that owns its physical state and uses Ballistic_new to compute pressure updates
        class Chamber {
        public:
            Chamber()
            {
                DefineConstants();
                pressure = 4.0;
            }
        
            void DefineConstants()
            {
                gamma = 1.4;
                At = 0.1;
                R = 287.0;
                T0 = 300.0;
                volume = 1.0;
            }
        
            void Update(Set::Scalar mdot, Set::Scalar dt)
            {
                pressure = ballistic.Compute(mdot, pressure, dt, At, gamma, R, T0, volume);
            }
        
            // Public state variables
            Set::Scalar pressure = NAN;
            Set::Scalar At = NAN;
            Set::Scalar gamma = NAN;
            Set::Scalar R = NAN;
            Set::Scalar T0 = NAN;
            Set::Scalar volume = NAN;
            Set::Scalar mdot = NAN;
        
        private:
            Ballistic_new ballistic;
        };

} // namespace Ballistics
} // namespace Model

#endif // MODEL_BALLISTICS_BALLISTIC_NEW_H_

#ifndef MODEL_INTERFACE_GB_H
#define MODEL_INTERFACE_GB_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace GB
{
class GB
{
    public:
    GB() {};
    virtual amrex::Real W(amrex::Real theta) = 0;
    virtual amrex::Real DW(amrex::Real theta) = 0;
    virtual amrex::Real DDW(amrex::Real theta) = 0;

    void ExportToFile(std::string filename, amrex::Real dTheta)
    {
        std::ofstream outFile;
        outFile.open(filename);
    
        for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
            {
                outFile << theta << " " << W(theta) << std::endl;   
            }
        outFile.close();
    
    }

protected:
    static constexpr amrex::Real pi = 3.14159265359;
};
}
}
}

#endif

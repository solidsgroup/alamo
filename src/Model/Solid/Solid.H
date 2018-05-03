#ifndef MODEL_SOLID_H_
#define MODEL_SOLID_H_

#include <AMReX.H>

namespace Model
{
namespace Solid
{

class Solid
{
public:
  Solid() {} ;

  virtual void 
  W(amrex::Real &w,
    const std::array<std::array<amrex::Real,2>,2> &gradu) = 0;

  virtual void
  DW(std::array<std::array<amrex::Real,2>,2> &dw,
     const std::array<std::array<amrex::Real,2>,2> &gradu) = 0;
  
  virtual void
  DDW(std::array<std::array<std::array<std::array<amrex::Real,2>,2>,2>,2> &ddw,
      const std::array<std::array<amrex::Real,2>,2> &gradu) = 0;

};

}
}


#endif



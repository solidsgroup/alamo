#ifndef MODEL_SOLID_H_
#define MODEL_SOLID_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{
class Solid
{
public:
  Solid() {} ;

  virtual Set::Scalar W(const Set::Matrix &gradu) = 0;

  virtual Set::Matrix DW(const Set::Matrix &gradu) = 0;
  virtual Set::Vector DW(const Set::Matrix &gradu, const Set::Vector &arg1) = 0;
  virtual Set::Scalar DW(const Set::Matrix &gradu, const Set::Matrix &arg1) = 0;

  //virtual Matrix DDW(const Matrix &gradu) = 0;
  virtual Set::Matrix DDW(const Set::Matrix &gradu, const Set::Vector &arg1, const Set::Vector &arg2) = 0;
  virtual Set::Scalar DDW(const Set::Matrix &gradu, const Set::Matrix &arg1, const Set::Matrix &arg2) = 0;

  // virtual void
  // DDW(std::array<std::array<std::array<std::array<amrex::Real,2>,2>,2>,2> &ddw,
  //     const std::array<std::array<amrex::Real,2>,2> &gradu) = 0;

};

}
}


#endif



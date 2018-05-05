#ifndef MODEL_SOLID_ELASTIC_H_
#define MODEL_SOLID_ELASTIC_H_

#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{

class Elastic : public Solid
{
public:
  Elastic();

  virtual Set::Scalar W(const Set::Matrix &gradu);

  virtual Set::Matrix DW(const Set::Matrix &gradu);
  virtual Set::Vector DW(const Set::Matrix &gradu, const Set::Vector &arg1);
  virtual Set::Scalar DW(const Set::Matrix &gradu, const Set::Matrix &arg1);

  //virtual Matrix DDW(const Matrix &gradu) = 0;
  virtual Set::Matrix DDW(const Set::Matrix &gradu, const Set::Vector &arg1, const Set::Vector &arg2);
  virtual Set::Scalar DDW(const Set::Matrix &gradu, const Set::Matrix &arg1, const Set::Matrix &arg2);

private:

  Set::Scalar mu, lambda;
};

}
}


#endif



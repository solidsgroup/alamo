#include "Elastic.H"

namespace Model
{
namespace Solid
{

Elastic::Elastic()
{
  mu = 1.0; lambda=1.0;
}

Set::Scalar Elastic::W(const  Set::Matrix &gradu)
{
  Set::Scalar w = 0.0;
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      w += mu * gradu(i,j) * gradu(k,l);
	    if (i==l && j==k) 
	      w += mu * gradu(i,j) * gradu(k,l);
	    if (i==j && k==l) 
	      w += lambda * gradu(i,j) * gradu(k,l);
	  }

  return 0.5*w;
}

Set::Matrix Elastic::DW(const Set::Matrix &gradu)
{
  Set::Matrix dw = Set::Matrix::Zero();
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      dw(i,j) += mu * gradu(k,l);
	    if (i==l && j==k) 
	      dw(i,j) += mu * gradu(k,l);
	    if (i==j && k==l) 
	      dw(i,j) += lambda * gradu(k,l);
	  }
  return dw;
}
Set::Vector Elastic::DW(const Set::Matrix &gradu, const Set::Vector &arg1)
{
  Set::Vector dw = Set::Vector::Zero();
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      dw(i) += mu * gradu(k,l)*arg1(j);
	    if (i==l && j==k) 
	      dw(i) += mu * gradu(k,l)*arg1(j);
	    if (i==j && k==l) 
	      dw(i) += lambda * gradu(k,l)*arg1(j);
	  }
  return dw;
}

Set::Scalar Elastic::DW(const Set::Matrix &gradu, const Set::Matrix &arg1)
{
  Set::Scalar dw = 0.0;
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      dw += mu * gradu(k,l)*arg1(i,j);
	    if (i==l && j==k) 
	      dw += mu * gradu(k,l)*arg1(i,j);
	    if (i==j && k==l) 
	      dw += lambda * gradu(k,l)*arg1(i,j);
	  }
  return dw;
}

Set::Matrix Elastic::DDW(const Set::Matrix & /*gradu*/, const Set::Vector &arg1, const Set::Vector &arg2)
{
  Set::Matrix ddw = Set::Matrix::Zero();
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      ddw(i,k) += mu * arg1(j) * arg2(l);
	    if (i==l && j==k) 
	      ddw(i,k) += mu * arg1(j) * arg2(l);
	    if (i==j && k==l) 
	      ddw(i,k) += lambda * arg1(j) * arg2(l);
	  }
  return ddw;
}

Set::Scalar Elastic::DDW(const Set::Matrix &/*gradu*/, const Set::Matrix &arg1, const Set::Matrix &arg2)
{
  Set::Scalar ddw = 0.0;
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<2; l++)
	  {
	    if (i==k && j==l) 
	      ddw += mu * arg1(i,j) * arg2(k,l);
	    if (i==l && j==k) 
	      ddw += mu * arg1(i,j) * arg2(k,l);
	    if (i==j && k==l) 
	      ddw += lambda * arg1(i,j) * arg2(k,l);
	  }
  return ddw;
}

}
}

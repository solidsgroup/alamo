#ifndef MODEL_INTERFACE_GB_COMBINE_H
#define MODEL_INTERFACE_GB_COMBINE_H

#include <iostream>
#include <fstream>
#include <vector>

#include "AMReX.H"
#include "GB.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
template<class GBE>
class Combine : public GB
{
public:
    Combine(int a_ngrains) : m_ngrains(a_ngrains) 
    {}
    Combine(IO::ParmParse &pp, int a_ngrains) : m_ngrains(a_ngrains) 
    {pp.queryclass(*this);}
    Combine(IO::ParmParse &pp, std::string name, int a_ngrains) : m_ngrains(a_ngrains) 
    {
        Util::Message(INFO);
        pp.queryclass(name,*this);
    }

    virtual Set::Scalar W(const Set::Scalar theta) const override
    {
        return m_w[0].W(theta);
        //Set::Scalar w = 0.0;
        //for (long unsigned int i = 0; i < m_mix.size(); i++) w += m_w[i].W(theta) * m_mix[i];
        //return w;
    };
    virtual Set::Scalar DW(const Set::Scalar theta) const
    {
        return m_w[0].DW(theta);
        //Set::Scalar dw = 0.0;
        //for (long unsigned int i = 0; i < m_mix.size(); i++) dw += m_w[i].DW(theta) * m_mix[i];
        //return dw;
    };
    virtual Set::Scalar DDW(const Set::Scalar theta) const
    {
        return m_w[0].DDW(theta);
        //Set::Scalar ddw = 0.0;
        //for (long unsigned int i = 0; i < m_mix.size(); i++) ddw += m_w[i].DDW(theta) * m_mix[i];
        //return ddw;
    };

    int m_ngrains = 0;
    std::vector<GBE> m_w;
    std::vector<Set::Scalar> m_mix;
    static constexpr Set::Scalar m_small = 1E-8;

    AMREX_FORCE_INLINE
    virtual void SetMix(amrex::Array4<const Set::Scalar> const &eta,int i, int j, int k) override
    {

        // Something going very wrong here. So we will abandon this for right now.
        //Set::Scalar etasqsum = 0.0;
        //int N = m_ngrains;
        //for (int I = 0 ; I < m_ngrains ; I++)
        //    for (int J = I+1 ; J < m_ngrains ; J++)
        //    {
        //        int idx = N*(N-1)/2 - (N-i)*(N-i-1)/2 + j - i - 1;
        //        m_mix[idx] = 1.0;//eta(i,j,k,I)*eta(i,j,k,I)*eta(i,j,k,J)*eta(i,j,k,J);
        //        etasqsum += m_mix[idx];
        //    }
        //for (int I = 0 ; I < m_mix.size(); I++)
        //{
        //    m_mix[I] /= (etasqsum + m_small);
        //}
        
    };

    static void Parse(Combine<GBE> & value, IO::ParmParse & pp)
    {
        value.m_w  .resize((value.m_ngrains) * (value.m_ngrains-1) / 2);
        value.m_mix.resize((value.m_ngrains) * (value.m_ngrains-1) / 2, NAN);

        Util::Message(INFO,value.m_ngrains);
        Util::Message(INFO,value.m_w.size());
        int N = value.m_ngrains;
        for (int i = 0 ; i < value.m_ngrains ; i++)
        for (int j = i+1 ; j < value.m_ngrains ; j++)
        {
            int idx = N*(N-1)/2 - (N-i)*(N-i-1)/2 + j - i - 1;
            Util::Message(INFO,"ngrains=",value.m_ngrains," i=",i," j=",j," idx=",idx);

            std::string name    = std::to_string(i)+"."+std::to_string(j);
            pp.queryclass(name.data(),value.m_w[idx]);            
        }
    }
      
};
}
}
}
#endif

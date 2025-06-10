#ifndef MODEL_REGRESSION_REGRESSION_H
#define MODEL_REGRESSION_REGRESSION_H

#include "Set/Base.H"
#include "IO/ParmParse.H"
#include <iterator>
#include <typeinfo>


namespace Model
{
namespace Regression
{

template<typename... OBJ>
class Regression
{
    friend IO::ParmParse;
private:
    std::tuple<OBJ...> obj;

    static constexpr const char* names[] = { OBJ::name... };
    
    int selected = -1;

public:

    template<int I = 0>
    void set_pressure (Set::Scalar P)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).set_pressure(P);
            else return set_pressure<I+1>(P);
        }
        else Util::Abort(INFO,selected);
    }
    
    template<int I = 0>
    Set::Scalar operator() (Set::Scalar phi, Set::Scalar T)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj)(phi,T);
            else return operator()<I+1>(phi,T);
        }
        else Util::Abort(INFO,selected);
    }


    template<int I = 0>
    static void Parse (Regression<OBJ...> &value, IO::ParmParse &pp)
    {
        if constexpr (I == 0)
        {
            std::string type;
            pp.query_default("type",type,value.names[0]);

            for (int i = 0; i < sizeof...(OBJ); i++)
                if (type == value.names[i])
                    value.selected = i;
            
            pp.pushPrefix(type);
        }

        if constexpr (I < sizeof...(OBJ))
        {
            if (value.selected == I)
            {
                std::get<I>(value.obj).Parse(std::get<I>(value.obj), pp);
                pp.popPrefix();
                return;
            }
            else 
                return Parse<I+1>(value,pp);
        }
        else Util::Abort(INFO);
    }
};






}
}




#endif

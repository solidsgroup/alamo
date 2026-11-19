#ifndef _TEST_SET_MATRIX4_H
#define _TEST_SET_MATRIX4_H

#include "Util/String.H"
#include "Util/Util.H"
#include "Set/Set.H"
#include "Set/Matrix4.H"

namespace Test
{
namespace Set
{
template <int dim, ::Set::Sym sym>
class Matrix4
{
public:
    static int PODTest(int verbose = 1)
    {
        int ret = 0;
        if (!std::is_trivially_copyable<::Set::Matrix4<dim,sym>>::value)
        {
            if (verbose) Util::Message(INFO,"Is not trivially copyable");
            ret += 1;
        }
        if (std::has_virtual_destructor<::Set::Matrix4<dim,sym>>::value)
        {
            if (verbose) Util::Message(INFO,"Should not have virtual destructor");
            ret += 1;
        }
        if (std::is_polymorphic<::Set::Matrix4<dim,sym>>::value)
        {
            if (verbose) Util::Message(INFO,"Should not be polymorphic");
            ret += 1;
        }
        if (!std::is_standard_layout<::Set::Matrix4<dim,sym>>::value)
        {
            if (verbose) Util::Message(INFO,"Does not have standard layout.");
            ret += 1;
        }
        Util::Test::SubMessage("POD", ret);
        return ret;
    }


    static int SymmetryTest(int /*verbose*/ = 0)
    {
        int failed = 0;
        ::Set::Matrix4<dim,sym> matrix = ::Set::Matrix4<dim,sym>::Random();

        //Util::Message(INFO,matrix);
        //std::cout << matrix << std::endl;
        if (sym == ::Set::Sym::MajorMinor)
        {
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                        {
                            if (matrix(i,j,k,l) != matrix(k,l,i,j)) failed += 1; // breaks major symmetry
                            if (matrix(i,j,k,l) != matrix(i,j,l,k)) failed += 1; // breaks minor symmetry
                            if (matrix(i,j,k,l) != matrix(j,i,k,l)) failed += 1; // breaks minor symmetry
                        }
        }
        else if (sym == ::Set::Sym::Full)
        {
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                        {
                            if (matrix(i,j,k,l) != matrix(i,j,l,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(i,k,j,l)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(i,k,l,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(i,l,j,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(i,l,k,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,i,k,l)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,i,l,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,k,i,l)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,k,l,i)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,l,i,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(j,l,k,i)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,i,j,l)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,i,l,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,j,i,l)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,j,l,i)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,l,i,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(k,l,j,i)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,i,j,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,i,k,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,j,i,k)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,j,k,i)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,k,i,j)) failed += 1;
                            if (matrix(i,j,k,l) != matrix(l,k,j,i)) failed += 1;
                        }
        }
        else 
        {
            return failed;
        }


        Util::Test::SubMessage("Symmetry", failed);
        return failed;
    }

    static int Test(int verbose = 1)
    {
        std::string funcname = __PRETTY_FUNCTION__;
        auto pos = funcname.find('[');
        if (pos != std::string_view::npos) funcname = funcname.substr(0, pos); 

        Util::Test::Message(funcname);
        int failed = 0;

        failed += PODTest(verbose);
        failed += SymmetryTest(verbose);
        return failed;
    }

};
}
}

#endif

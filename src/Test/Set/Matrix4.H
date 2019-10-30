#include "Set/Set.H"
namespace Test
{
namespace Set
{
template <int dim, int sym>
class Matrix4
{
public:
    int SymmetryTest(int /*verbose*/)
    {
        ::Set::Matrix4<dim,sym> matrix = ::Set::Matrix4<dim,sym>::Randomize();

        //Util::Message(INFO,matrix);
        //std::cout << matrix << std::endl;
        if (sym == ::Set::Sym::None)
            return 0;
        if (sym == ::Set::Sym::MajorMinor)
        {
            for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            for (int k = 0; k < dim; k++)
            for (int l = 0; l < dim; l++)
            {
                if (matrix(i,j,k,l) != matrix(k,l,i,j)) return 1; // breaks major symmetry
                if (matrix(i,j,k,l) != matrix(i,j,l,k)) return 1; // breaks minor symmetry
                if (matrix(i,j,k,l) != matrix(j,i,k,l)) return 1; // breaks minor symmetry
            }
            return 0;
        }
        if (sym == ::Set::Sym::Full)
        {
            for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            for (int k = 0; k < dim; k++)
            for (int l = 0; l < dim; l++)
            {
                if (matrix(i,j,k,l) != matrix(i,j,l,k)) return 1;
                if (matrix(i,j,k,l) != matrix(i,k,j,l)) return 1;
                if (matrix(i,j,k,l) != matrix(i,k,l,j)) return 1;
                if (matrix(i,j,k,l) != matrix(i,l,j,k)) return 1;
                if (matrix(i,j,k,l) != matrix(i,l,k,j)) return 1;
                if (matrix(i,j,k,l) != matrix(j,i,k,l)) return 1;
                if (matrix(i,j,k,l) != matrix(j,i,l,k)) return 1;
                if (matrix(i,j,k,l) != matrix(j,k,i,l)) return 1;
                if (matrix(i,j,k,l) != matrix(j,k,l,i)) return 1;
                if (matrix(i,j,k,l) != matrix(j,l,i,k)) return 1;
                if (matrix(i,j,k,l) != matrix(j,l,k,i)) return 1;
                if (matrix(i,j,k,l) != matrix(k,i,j,l)) return 1;
                if (matrix(i,j,k,l) != matrix(k,i,l,j)) return 1;
                if (matrix(i,j,k,l) != matrix(k,j,i,l)) return 1;
                if (matrix(i,j,k,l) != matrix(k,j,l,i)) return 1;
                if (matrix(i,j,k,l) != matrix(k,l,i,j)) return 1;
                if (matrix(i,j,k,l) != matrix(k,l,j,i)) return 1;
                if (matrix(i,j,k,l) != matrix(l,i,j,k)) return 1;
                if (matrix(i,j,k,l) != matrix(l,i,k,j)) return 1;
                if (matrix(i,j,k,l) != matrix(l,j,i,k)) return 1;
                if (matrix(i,j,k,l) != matrix(l,j,k,i)) return 1;
                if (matrix(i,j,k,l) != matrix(l,k,i,j)) return 1;
                if (matrix(i,j,k,l) != matrix(l,k,j,i)) return 1;
            }
            return 0;
        }


        return 1;
    }
};
}
}
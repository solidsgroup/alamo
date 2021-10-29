
inline void operator+= (amrex::Vector<amrex::Real> a)
{
    eigen_assert(a.size() == AMREX_SPACEDIM);
    AMREX_D_TERM(this->operator()(0) += a[0];,
            this->operator()(1) += a[1];,
            this->operator()(2) += a[2];);
}
Matrix<amrex::Real,AMREX_SPACEDIM,1> operator* (const amrex::Vector<amrex::Real>& x) // WORKS
{
    eigen_assert(x.size() == AMREX_SPACEDIM);
    Matrix<amrex::Real,AMREX_SPACEDIM,1> ret;
    AMREX_D_TERM(ret(0) = this->operator()(0) * x[0];,
            ret(1) = this->operator()(1) * x[1];,
            ret(2) = this->operator()(2) * x[2];);
    return ret;
}
Matrix<_Scalar,_Rows,_Cols>& operator*= (const amrex::Vector<amrex::Real>& x) // WORKS
{
    eigen_assert(x.size() == AMREX_SPACEDIM);
    AMREX_D_TERM(this->operator()(0) *= x[0];,
            this->operator()(1) *= x[1];,
            this->operator()(2) *= x[2];);
    return *this;
}
operator amrex::Vector<amrex::Real>() // WORKS
{
    return {{AMREX_D_DECL(this->operator()(0),
                this->operator()(1),
                this->operator()(2))}};
}

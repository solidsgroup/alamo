AMREX_FORCE_INLINE
void operator += (const OP_CLASS &rhs)
{
#define X(name) \
    name += rhs.name;
    OP_VARS
#undef X
}    

AMREX_FORCE_INLINE
OP_CLASS operator * (const Set::Scalar alpha) const
{
        OP_CLASS ret;
#define X(name) \
    ret.name = name*alpha;
    OP_VARS
#undef X
    return ret;
}    
friend OP_CLASS operator * (const Set::Scalar alpha, const OP_CLASS b);
friend OP_CLASS operator + (const OP_CLASS a, const OP_CLASS b);
friend OP_CLASS operator - (const OP_CLASS a, const OP_CLASS b);

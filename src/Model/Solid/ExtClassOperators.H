AMREX_FORCE_INLINE
OP_CLASS operator * (const Set::Scalar alpha, const OP_CLASS b)
{
    OP_CLASS ret;
#define X(name) \
    ret.name = b.name*alpha;
    OP_VARS
#undef X
    return ret;
}    

AMREX_FORCE_INLINE
OP_CLASS operator + (const OP_CLASS a, const OP_CLASS b)
{
    OP_CLASS ret;
#define X(name) \
    ret.name = a.name + b.name;
    OP_VARS
#undef X
    return ret;
}    

AMREX_FORCE_INLINE
OP_CLASS operator - (const OP_CLASS a, const OP_CLASS b)
{
    OP_CLASS ret;
#define X(name) \
    ret.name = a.name - b.name;
    OP_VARS
#undef X
    return ret;
}    

AMREX_FORCE_INLINE
bool operator == (const OP_CLASS a, const OP_CLASS b)
{
    bool ret = true;
#define X(name) \
    if (!(a.name == b.name)) ret = false;
    OP_VARS
#undef X
    return ret;
}

#undef OP_VARS
#undef OP_CLASS

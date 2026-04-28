#ifndef UTIL_NB_UTIL_H_
#define UTIL_NB_UTIL_H_

namespace NarrowBandData
{
    struct TubeIdx
    {
        static constexpr int NB_MASK   = 0;  // 0/1: in narrow band
        static constexpr int ZERO_ID   = 1;  // zero level set id
        static constexpr int OBJ_ID    = 2;  // object id
        static constexpr int TUBE_TYPE = 3;  // enum below

        static constexpr int NCOMP     = 4;
    };

    struct TubeType
    {
        static constexpr int INTERFACE   = 0;
        static constexpr int INNERBAND   = 1;
        static constexpr int INSIDETUBE  = 2;
        static constexpr int EDGEPOINT   = 3;
        static constexpr int OUTSIDETUBE = 4;
    };
};

#endif
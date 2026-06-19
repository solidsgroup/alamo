#ifndef UTIL_BMP_H
#define UTIL_BMP_H

#include <limits>
#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "Util/Util.H"

namespace Util
{
class BMP
{
public:
    BMP () {};

    void Define (std::string filename)
    {
        FILE* f;
        int retcode = -1;
        if ((f = fopen(filename.c_str(), "rb")))
        {
            unsigned char info[54];
            retcode = fread(info, sizeof(unsigned char), 54, f); 
            if (retcode != 54) Util::Abort(INFO,"Error reading BMP file: 54 bytes expected but only ", retcode, " read.");
            nx = *(int*)&info[18];
            ny = *(int*)&info[22];
            
            std::vector<std::array<int,3>> data_host;
            data_host.resize(nx*ny);
            
            int row_padded = (nx*3 + 3) & (~3);
            unsigned char* chardata = new unsigned char[row_padded];
            for (int j = 0; j < ny; j++)
            {
                retcode = fread(chardata, sizeof(unsigned char), row_padded, f);
                if (retcode != row_padded) Util::Abort(INFO,"Error reading BMP file: ",row_padded," bytes expected but only ", retcode, " read.");
                for(int i = 0; i < nx*3; i += 3)
                {
                    data_host[i/3 + nx*j][0] = (int)chardata[i+2]; // R
                    data_host[i/3 + nx*j][1] = (int)chardata[i+1]; // G
                    data_host[i/3 + nx*j][2] = (int)chardata[i];   // B

                    m_min[0] = std::min(m_min[0], data_host[i/3 + nx*j][0]);
                    m_min[1] = std::min(m_min[1], data_host[i/3 + nx*j][1]);
                    m_min[2] = std::min(m_min[2], data_host[i/3 + nx*j][2]);

                    m_max[0] = std::min(m_max[0], data_host[i/3 + nx*j][0]);
                    m_max[1] = std::min(m_max[1], data_host[i/3 + nx*j][1]);
                    m_max[2] = std::min(m_max[2], data_host[i/3 + nx*j][2]);
                }
            }
            fclose(f);
            delete[] chardata;

            data.resize(data_host.size());
            amrex::Gpu::copyAsync(  amrex::Gpu::hostToDevice,
                                    data_host.begin(),
                                    data_host.end(),
                                    data.begin());
            amrex::Gpu::streamSynchronize();
            data_ptr = data.data();
        }
        else
        {
            Util::Abort(INFO,"File ", filename, " does not exist");
        }
    }

    AMREX_FORCE_INLINE
    AMREX_GPU_DEVICE
    const std::array<int,3> operator () (int i,int j) const
    {
        Util::Assert(INFO,TEST(i < nx)," i = ",i," nx = ", nx);
        Util::Assert(INFO,TEST(j < ny)," j = ",j," ny = ", ny);
        return data_ptr[nx*j + i];
    }

    AMREX_GPU_HOST
    std::array<int,3> min()
    {
        return m_min;
    }
    AMREX_GPU_HOST
    std::array<int,3> max()
    {
        return m_max;
    }

public:
    int nx, ny;
    std::array<int, 3> m_min = {
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max()
    };
    std::array<int, 3> m_max = {
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min()};

    amrex::Gpu::DeviceVector<std::array<int,3>> data;
    std::array<int,3> *data_ptr = nullptr;
};
}


#endif

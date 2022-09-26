#ifndef UTIL_BMP_H
#define UTIL_BMP_H

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
            
            data.resize(nx*ny);
            
            int row_padded = (nx*3 + 3) & (~3);
            unsigned char* data = new unsigned char[row_padded];
            for (int j = 0; j < ny; j++)
            {
                retcode = fread(data, sizeof(unsigned char), row_padded, f);
                if (retcode != row_padded) Util::Abort(INFO,"Error reading BMP file: ",row_padded," bytes expected but only ", retcode, " read.");
                for(int i = 0; i < nx*3; i += 3)
                {
                    (*this)(i/3,j)[0] = (int)data[i+2]; // R
                    (*this)(i/3,j)[1] = (int)data[i+1]; // G
                    (*this)(i/3,j)[2] = (int)data[i];   // B
                }
            }
            fclose(f); 
        }
        else
        {
            Util::Abort(INFO,"File ", filename, " does not exist");
        }
    }

    inline
    std::array<int,3> & operator () (int i,int j)
    {
        Util::Assert(INFO,TEST(i < nx)," i = ",i," nx = ", nx);
        Util::Assert(INFO,TEST(j < ny)," j = ",j," ny = ", ny);
        return data[nx*j + i];
    }
    std::array<int,3> min()
    {
        std::array<int,3> _min = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        for (unsigned int i = 0; i < data.size(); i++)
        {
            if (_min[0] > data[i][0]) _min[0] = data[i][0];
            if (_min[1] > data[i][1]) _min[1] = data[i][1];
            if (_min[2] > data[i][2]) _min[2] = data[i][2];
        }
        return _min;
    }
    std::array<int,3> max()
    {
        std::array<int,3> _max = {0, 0, 0};
        for (unsigned int i = 0; i < data.size(); i++)
        {
            if (_max[0] < data[i][0]) _max[0] = data[i][0];
            if (_max[1] < data[i][1]) _max[1] = data[i][1];
            if (_max[2] < data[i][2]) _max[2] = data[i][2];
        }
        return _max;
    }

public:
    int nx, ny;
private:
    std::vector<std::array<int,3>> data;
};
}


#endif

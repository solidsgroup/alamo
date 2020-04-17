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
        FILE* f = fopen(filename.c_str(), "rb");
        unsigned char info[54];
        fread(info, sizeof(unsigned char), 54, f); 
        nx = *(int*)&info[18];
        ny = *(int*)&info[22];
        Util::Message(INFO);
        
        data.resize(nx*ny);
        
        int row_padded = (nx*3 + 3) & (~3);
        unsigned char* data = new unsigned char[row_padded];
        for (int i = 0; i < ny; i++)
        {
            fread(data, sizeof(unsigned char), row_padded, f);
            for(int j = 0; j < nx*3; j += 3)
            {
                (*this)(j/3,i)[0] = (int)data[j+2]; // R
                (*this)(j/3,i)[1] = (int)data[j+1]; // G
                (*this)(j/3,i)[2] = (int)data[j];   // B
            }
        }
        fclose(f); 
    }

    inline
    std::array<int,3> & operator () (int i,int j)
    {
        Util::Assert(INFO,TEST(i < nx),"i = ",i," nx = ", nx);
        Util::Assert(INFO,TEST(j < ny),"j = ",j," ny = ", ny);
        return data[nx*j + i];
    }

public:
    int nx, ny;
private:
    std::vector<std::array<int,3>> data;
};
}


#endif
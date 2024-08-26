#ifndef UTIL_BMP4_H  // Updated include guard to match the new filename
#define UTIL_BMP4_H

#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <limits>
#include "Util/Util.H"

namespace Util
{
class BMP4
{
public:
    BMP4() : nx(0), ny(0) {}  // Initialize nx and ny to 0

    void Define(std::string filename)
    {
        FILE* f;
        int retcode = -1;
        if ((f = fopen(filename.c_str(), "rb")))
        {
            unsigned char info[54];
            retcode = fread(info, sizeof(unsigned char), 54, f); 
            if (retcode != 54) 
                Util::Abort(INFO, "Error reading BMP file: 54 bytes expected but only ", retcode, " read.");
            
            nx = *(int*)&info[18];  // Width of the BMP image
            ny = *(int*)&info[22];  // Height of the BMP image
            
            data.resize(nx * ny);  // Resize the data vector to hold nx*ny pixels
            
            int row_padded = (nx * 4 + 3) & (~3); // Adjusted for 4 channels (RGBA), ensures row alignment to 4 bytes
            unsigned char* buffer = new unsigned char[row_padded]; // Buffer to hold one row of pixel data
            
            for (int j = 0; j < ny; j++)
            {
                retcode = fread(buffer, sizeof(unsigned char), row_padded, f);
                if (retcode != row_padded) 
                    Util::Abort(INFO, "Error reading BMP file: ", row_padded, " bytes expected but only ", retcode, " read.");
                
                for (int i = 0; i < nx * 4; i += 4)  // Adjusted for 4 channels (RGBA)
                {
                    // Assign each channel (R, G, B, A) from the buffer to the appropriate place in the data array
                    (*this)(i / 4, j)[0] = (int)buffer[i + 2]; // R
                    (*this)(i / 4, j)[1] = (int)buffer[i + 1]; // G
                    (*this)(i / 4, j)[2] = (int)buffer[i];     // B
                    (*this)(i / 4, j)[3] = (int)buffer[i+3]; // A (Alpha channel)
                }
            }
            delete[] buffer; // Clean up the allocated buffer memory
            fclose(f);
        }
        else
        {
            Util::Abort(INFO, "File ", filename, " does not exist");
        }
    }

    AMREX_FORCE_INLINE
    std::array<int, 4>& operator()(int i, int j)  // Now returns a reference to an array of 4 integers
    {
        Util::Assert(INFO, TEST(i < nx), " i = ", i, " nx = ", nx);
        Util::Assert(INFO, TEST(j < ny), " j = ", j, " ny = ", ny);
        return data[nx * j + i];  // Access the correct pixel (i, j) in the flattened data vector
    }

    std::array<int, 4> min()  // Finds the minimum value for each of the 4 channels
    {
        std::array<int, 4> _min = {
            std::numeric_limits<int>::max(), 
            std::numeric_limits<int>::max(), 
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()  // Initialize with the maximum possible int value
        };

        for (unsigned int i = 0; i < data.size(); i++)
        {
            for (int c = 0; c < 4; ++c)  // Iterate over each channel (R, G, B, A)
                if (_min[c] > data[i][c]) _min[c] = data[i][c];  // Update the minimum value for each channel
        }
        return _min;
    }

    std::array<int, 4> max()  // Finds the maximum value for each of the 4 channels
    {
        std::array<int, 4> _max = {0, 0, 0, 0};  // Initialize with zeroes for all 4 channels
        for (unsigned int i = 0; i < data.size(); i++)
        {
            for (int c = 0; c < 4; ++c)  // Iterate over each channel (R, G, B, A)
                if (_max[c] < data[i][c]) _max[c] = data[i][c];  // Update the maximum value for each channel
        }
        return _max;
    }

public:
    int nx, ny;  // Dimensions of the BMP image

private:
    std::vector<std::array<int, 4>> data;  // Now stores an array of 4 integers per pixel (R, G, B, A)
};
}

#endif

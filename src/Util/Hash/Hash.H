#ifndef UTIL_HASH_HASH_H
#define UTIL_HASH_HASH_H

#include <iostream>
#include <fstream>
#include <regex>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>

namespace Util
{
namespace Hash
{

// Helper function to check if a string ends with a specific suffix
bool endsWith(const std::string& str, const std::string& suffix) {
    if (str.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

// Function to remove comments and extra whitespace from a file
std::string preprocessSourceFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file) {
        std::cerr << "Cannot open file: " << filePath << std::endl;
        return "";
    }

    std::ostringstream preprocessed;
    std::string line;
    bool inBlockComment = false;

    while (std::getline(file, line)) {
        // Remove block comments
        if (inBlockComment) {
            size_t endBlock = line.find("*/");
            if (endBlock != std::string::npos) {
                line = line.substr(endBlock + 2);
                inBlockComment = false;
            } else {
                continue;  // Skip lines inside block comments
            }
        }

        // Remove line comments
        size_t lineComment = line.find("//");
        if (lineComment != std::string::npos) {
            line = line.substr(0, lineComment);
        }

        // Remove block comments that start on the same line
        size_t startBlock = line.find("/*");
        if (startBlock != std::string::npos) {
            size_t endBlock = line.find("*/", startBlock + 2);
            if (endBlock != std::string::npos) {
                line.erase(startBlock, endBlock - startBlock + 2);
            } else {
                line.erase(startBlock);
                inBlockComment = true;
            }
        }

        // Trim extra spaces (optional)
        std::regex extraSpace("\\s+");
        line = std::regex_replace(line, extraSpace, " ");

        // Append the cleaned line if it’s not empty
        if (!line.empty()) {
            preprocessed << line << "\n";
        }
    }

    return preprocessed.str();
}

// Native implementation of SHA-256
class SHA256 {
public:
    SHA256() { reset(); }

    void update(const uint8_t* data, size_t length) {
        for (size_t i = 0; i < length; ++i) {
            data_[datalen_++] = data[i];
            if (datalen_ == 64) {
                transform();
                bitlen_ += 512;
                datalen_ = 0;
            }
        }
    }

    void update(const std::string& data) {
        update(reinterpret_cast<const uint8_t*>(data.c_str()), data.size());
    }

    std::string final() {
        uint8_t hash[32];
        pad();
        transform();
        toBytes(hash, state_, 8);

        std::ostringstream result;
        for (int i = 0; i < 32; ++i) {
            result << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
        }

        return result.str();
    }

    void reset() {
        datalen_ = 0;
        bitlen_ = 0;
        state_ = {  0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
    }

private:
    std::array<uint32_t, 8> state_;
    uint8_t data_[64] = {0};
    uint32_t datalen_ = 0;
    uint64_t bitlen_ = 0;

    const std::array<uint32_t, 64> k_ = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b,
        0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01,
        0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7,
        0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152,
        0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
        0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819,
        0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08,
        0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f,
        0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

    static void toBytes(uint8_t* output, const std::array<uint32_t, 8>& input, size_t size) {
        for (size_t i = 0; i < size; ++i) {
            output[i * 4 + 0] = (input[i] >> 24) & 0xff;
            output[i * 4 + 1] = (input[i] >> 16) & 0xff;
            output[i * 4 + 2] = (input[i] >> 8) & 0xff;
            output[i * 4 + 3] = (input[i] >> 0) & 0xff;
        }
    }

    static uint32_t rotr(uint32_t x, uint32_t n) {
        return (x >> n) | (x << (32 - n));
    }

    static uint32_t choose(uint32_t e, uint32_t f, uint32_t g) {
        return (e & f) ^ (~e & g);
    }

    static uint32_t majority(uint32_t a, uint32_t b, uint32_t c) {
        return (a & b) ^ (a & c) ^ (b & c);
    }

    static uint32_t sig0(uint32_t x) {
        return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3);
    }

    static uint32_t sig1(uint32_t x) {
        return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10);
    }

    void transform() {
        std::array<uint32_t, 64> m;
        std::array<uint32_t, 8> v = state_;

        for (size_t i = 0, j = 0; i < 16; ++i, j += 4) {
            m[i] = (data_[j] << 24) | (data_[j + 1] << 16) | (data_[j + 2] << 8) | data_[j + 3];
        }

        for (size_t i = 16; i < 64; ++i) {
            m[i] = sig1(m[i - 2]) + m[i - 7] + sig0(m[i - 15]) + m[i - 16];
        }

        for (size_t i = 0; i < 64; ++i) {
            uint32_t temp1 = v[7] + choose(v[4], v[5], v[6]) + rotr(v[4], 6) + k_[i] + m[i];
            uint32_t temp2 = majority(v[0], v[1], v[2]) + rotr(v[0], 2);
            v[7] = v[6];
            v[6] = v[5];
            v[5] = v[4];
            v[4] = v[3] + temp1;
            v[3] = v[2];
            v[2] = v[1];
            v[1] = v[0];
            v[0] = temp1 + temp2;
        }

        for (size_t i = 0; i < 8; ++i) {
            state_[i] += v[i];
        }
    }

    void pad() {
        size_t i = datalen_;

        if (datalen_ < 56) {
            data_[i++] = 0x80;
            while (i < 56) {
                data_[i++] = 0x00;
            }
        } else {
            data_[i++] = 0x80;
            while (i < 64) {
                data_[i++] = 0x00;
            }
            transform();
            std::fill(std::begin(data_), std::end(data_), 0);
        }

        bitlen_ += datalen_ * 8;
        data_[63] = bitlen_;
        data_[62] = bitlen_ >> 8;
        data_[61] = bitlen_ >> 16;
        data_[60] = bitlen_ >> 24;
        data_[59] = bitlen_ >> 32;
        data_[58] = bitlen_ >> 40;
        data_[57] = bitlen_ >> 48;
        data_[56] = bitlen_ >> 56;
    }
};

// Function to process and hash a single file
std::string processFile(const std::string& filePath) {
    std::string preprocessed = preprocessSourceFile(filePath);
    if (preprocessed.empty()) {
        return "";
    }

    SHA256 sha256;
    sha256.update(preprocessed);
    return sha256.final();
}

// Recursive function to scan all files in a directory, process .cpp and .H files, and sort them
void processDirectory(const std::string& directory, std::vector<std::string>& sortedFiles) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            std::string filePath = entry.path().string();
            if (endsWith(filePath, ".cpp") || endsWith(filePath, ".H")) {
                sortedFiles.push_back(filePath);
            }
        }
    }

    // Sort files in lexicographical order
    std::sort(sortedFiles.begin(), sortedFiles.end());
}

int getFinalHash(std::string srcDirectory) {
    std::vector<std::string> sortedFiles;

    std::cout << "Scanning directory: " << srcDirectory << std::endl;
    processDirectory(srcDirectory, sortedFiles);

    // Process and hash all files together
    SHA256 sha256;
    for (const auto& filePath : sortedFiles) {
        std::cout << "Processing file: " << filePath << std::endl;
        std::string fileHash = processFile(filePath);
        sha256.update(fileHash);
    }

    std::cout << "\nFinal Hash (SHA256 of all preprocessed content): " << sha256.final() << std::endl;

    return 0;
}



}
}

#endif

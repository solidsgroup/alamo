#ifndef UTIL_STRING_H
#define UTIL_STRING_H

#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>
#include <complex>
#include <stdlib.h>

#include <string>

#include "AMReX_Config.H"
#include "Set/Base.H"
#include "Util/Util.H"

/// \brief A collection of utility routines
namespace Util
{
/// \brief A collection of string operations
namespace String
{
AMREX_FORCE_INLINE
int ReplaceAll(std::string &str, const std::string before, const std::string after)
{
    size_t start_pos = 0;
    while((start_pos = str.find(before, start_pos)) != std::string::npos) {
        str.replace(start_pos, before.length(), after);
        start_pos += after.length();
    }
    return 0;
}

AMREX_FORCE_INLINE
int ReplaceAll(std::string &str, const char before, const std::string after)
{
    size_t start_pos = 0;
    while((start_pos = str.find(before, start_pos)) != std::string::npos) {
        str.replace(start_pos, 1, after);
        start_pos += after.length();
    }
    return 0;
}

AMREX_FORCE_INLINE
std::string Wrap(std::string text, unsigned per_line)
{
    unsigned line_begin = 0;

    while (line_begin < text.size())
    {
        const unsigned ideal_end = line_begin + per_line ;
        unsigned line_end = ideal_end <= text.size() ? ideal_end : text.size()-1;

        if (line_end == text.size() - 1)
            ++line_end;
        else if (std::isspace(text[line_end]))
        {
            text[line_end] = '\n';
            ++line_end;
        }
        else    // backtrack
        {
            unsigned end = line_end;
            while ( end > line_begin && !std::isspace(text[end]))
                --end;

            if (end != line_begin)                  
            {                                       
                line_end = end;                     
                text[line_end++] = '\n';            
            }                                       
            else                                    
                text.insert(line_end++, 1, '\n');
        }

        line_begin = line_end;
    }

    return text;
}

AMREX_FORCE_INLINE
std::string Join(const std::vector<std::string> & vec, std::string separator = "_")
{
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i < vec.size() - 1) { // Don't add separator after the last element
            oss << separator;
        }
    }
    return oss.str();
}

AMREX_FORCE_INLINE
std::string Join(const std::vector<std::string> & vec, char separator)
{
    return Join(vec,std::to_string(separator));
}

AMREX_FORCE_INLINE
std::string Join(const std::vector<double> & vec, std::string separator = " ")
{
    std::vector<std::string> strvec;
    for (auto v : vec) strvec.push_back(std::to_string(v));
    return Join(strvec,separator);
}

AMREX_FORCE_INLINE
std::string Join(const Set::Vector & vec, char separator = ' ')
{
    std::ostringstream oss;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        oss << vec(i);
        if (i < vec.size() - 1) { // Don't add separator after the last element
            oss << separator;
        }
    }
    return oss.str();
}

AMREX_FORCE_INLINE
std::string Join(const Set::Matrix & vec, char separator = ' ')
{
    std::ostringstream oss;
    for (unsigned int i = 0; i < AMREX_SPACEDIM; ++i) {
        for (unsigned int j = 0; j < AMREX_SPACEDIM; ++j) {
            oss << vec(i,j);
            if (i < AMREX_SPACEDIM - 1 && j < AMREX_SPACEDIM - 1 ) { // Don't add separator after the last element
                oss << separator;
            }
        }
    }
    return oss.str();
}

AMREX_FORCE_INLINE
std::vector<std::string> Split(std::string &str, const char delim = ' ')
{
    std::vector<std::string> ret;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        ret.push_back(item);
    }
    return ret;
}

AMREX_FORCE_INLINE
bool Contains(std::string &str, const std::string find)
{
    if (str.find(find) != std::string::npos) return true;
    else return false;
}

template<class T>
AMREX_FORCE_INLINE
T Parse(std::string);

template<>
AMREX_FORCE_INLINE
std::complex<int> Parse(std::string input)
{
    int re=0, im=0;

    ReplaceAll(input, "+", " +");
    ReplaceAll(input, "-", " -");
    std::vector<std::string> tokens = Split(input,' ');
    for (unsigned int i = 0; i < tokens.size(); i++)
    {
        if(tokens[i]=="") continue;
        if(Contains(tokens[i],"i"))
        {
            ReplaceAll(tokens[i],"i","");
            im += std::stoi(tokens[i]);
        }
        else 
        {
            re += std::stoi(tokens[i]);
        }
    }
    return std::complex<int>(re,im);
}

}

}
#endif

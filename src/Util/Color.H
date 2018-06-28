#ifndef UTIL_COLORS_H
#define UTIL_COLORS_H

#include <string>

namespace Color
{
static std::string Reset        = "\033[0m";
static std::string Bold         = "\033[1m";
static std::string Dim          = "\033[2m";
static std::string Underlined   = "\033[4m";
static std::string Blink        = "\033[5m";
static std::string Reverse      = "\033[7m";
static std::string Hiden        = "\033[8m";

namespace FG
{
static std::string Default      = "\033[39m";
static std::string Black        = "\033[30m";
static std::string Red          = "\033[31m";
static std::string Green        = "\033[32m";
static std::string Yellow       = "\033[33m";
static std::string Blue         = "\033[34m";
static std::string Magenta      = "\033[35m";
static std::string Cyan         = "\033[36m";
static std::string LightGray    = "\033[37m";
static std::string DarkGray     = "\033[90m";
static std::string LightRed     = "\033[91m";
static std::string LightGreen   = "\033[92m";
static std::string LightYellow  = "\033[93m";
static std::string LightBlue    = "\033[94m";
static std::string LightMagenta = "\033[95m";
static std::string LightCyan    = "\033[96m";
static std::string White        = "\033[97m";
}

namespace BG
{
static std::string Default      = "\033[49m";
static std::string Red          = "\033[41m";
static std::string Green        = "\033[42m";
static std::string Yellow       = "\033[43m";
static std::string Blue         = "\033[44m";
static std::string Magenta      = "\033[45m";
static std::string Cyan         = "\033[46m";
static std::string LightGray    = "\033[47m";
static std::string DarkGray     = "\033[100m";
static std::string LightRed     = "\033[101m";
static std::string LightGreen   = "\033[102m";
static std::string LightYellow  = "\033[103m";
static std::string LightBlue    = "\033[104m";
static std::string LightMagenta = "\033[105m";
static std::string LightCyan    = "\033[106m";
static std::string White        = "\033[107m";
}

}

#endif//WIELD_UTILS_COLOR_H

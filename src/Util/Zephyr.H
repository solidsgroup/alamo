#ifndef UTIL_ZEPHYR_H
#define UTIL_ZEPHYR_H

#include <string>

namespace Util::Zephyr
{
void Enable();
void Start(const std::string& output_directory);
void Stop();
}

#endif

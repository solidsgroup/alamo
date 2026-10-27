#ifndef IO_OUTPUTLOG_H
#define IO_OUTPUTLOG_H

#include <string>

namespace IO::OutputLog
{

/// Route stdout and stderr through the terminal/file tee.
void Initialize();

/// Open the rank-zero log after the output directory has been created.
void Open(const std::string& path, bool append = false);

/// Stop retaining output on ranks that do not own the log file.
void DisableFile();

/// Flush the log and restore the original terminal stream buffers.
void Finalize();

}

#endif

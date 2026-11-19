//
// ============================================
// :fas:`file-code;fa-fw` Filename parsing
// ============================================
//
// The :code:`plot_file` input defines the simulation output directory. Alamo
// supports time and dimension wildcards, input-variable substitution, and
// Slurm environment substitution in this value.
//
// Time and dimension wildcards
// ----------------------------
//
// The following wildcard strings are available:
//
// .. code-block:: text
//
//    %Y  Four-digit year
//    %m  Two-digit month (01-12)
//    %d  Two-digit day of the month (01-31)
//    %H  Two-digit hour (00-23)
//    %M  Two-digit minute (00-59)
//    %S  Two-digit second (00-59)
//    %f  Microseconds within the current second
//    %D  Spatial dimension
//
// For example:
//
// .. code-block:: makefile
//
//    plot_file = output_%Y%m%d_%H%M%S
//
// may resolve to:
//
// .. code-block:: text
//
//    output_20250307_211201
//
// Input variable substitution
// ---------------------------
//
// Use braces to substitute a variable from the input file:
//
// .. code-block:: makefile
//
//    myinput.value = 3
//    plot_file = output_{myinput.value}
//
// This resolves the output directory to:
//
// .. code-block:: text
//
//    output_3
//
// Slurm environment substitution
// ------------------------------
//
// Any brace-delimited variable whose name begins with :code:`SLURM_` is read
// from the process environment. Common substitutions include:
//
// .. code-block:: text
//
//    {SLURM_JOB_ID}
//    {SLURM_JOB_NAME}
//    {SLURM_JOB_PARTITION}
//    {SLURM_CLUSTER_NAME}
//    {SLURM_ARRAY_TASK_ID}
//
// For example:
//
// .. code-block:: makefile
//
//    plot_file = output_{SLURM_JOB_ID}
//
// may resolve to:
//
// .. code-block:: text
//
//    output_314159
//
// If the requested variable is unavailable or empty, its name is retained
// without braces:
//
// .. code-block:: text
//
//    output_SLURM_JOB_ID

#ifndef IO_FILENAMEPARSE_H
#define IO_FILENAMEPARSE_H

#include <string>


namespace IO
{

void FileNameParse(std::string &filename);

}

#endif

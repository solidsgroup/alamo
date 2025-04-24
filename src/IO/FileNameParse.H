//
// This manages the processing of the :code:`plot_file` input, which defines the output file.
// Both variable substitution and time wildcards are supported.
//
// **Variable substitution**: you can use braces to substitute variable names in the input file.
// For example, if your input file contains :code:`myinput.value = 3`, then you can specify
// :code:`plot_file = output_{myinput.value}` to substitute the input's value.
//
// **Time and dimension wildcards**: You can use the following strings:
//
// - :code:`%Y`: 4-digit year
// - :code:`%m`: 2-digit month (01-12)
// - :code:`%d`: 2-digit day (00-31)
// - :code:`%H`: 2 digit hour (00-23)
// - :code:`%M`: 2 digit minute (00-59)
// - :code:`%S`: 2 digit second (00-59)
// - :code:`%f`: 2 digit microsecond (00-59)
// - :code:`%D`: spatial dimension
//
// So a :code:`plot_file = output_%Y%m%d_%H%M%S` may resolve to
// :code:`output_20250307_211201` (etc) depending on the time the
// output was created.

#ifndef IO_FILENAMEPARSE_H
#define IO_FILENAMEPARSE_H

#include <string>


namespace IO
{

void FileNameParse(std::string &filename);

}

#endif

#include "ParmParse.H"

namespace IO
{
bool ParmParse::checked_for_input_files = false;

void
ParmParse::SetTraversalMode(bool enabled)
{
    InputScraper::SetTraversalMode(enabled);
}

bool
ParmParse::InTraversalMode()
{
    return InputScraper::InTraversalMode();
}

bool
ParmParse::ShouldExecute()
{
    return InputScraper::ShouldExecute();
}

const ParmParse::InputNode &
ParmParse::InputTree()
{
    return InputScraper::InputTree();
}

void
ParmParse::ClearInputTree()
{
    InputScraper::ClearInputTree();
}

void
ParmParse::SetTraversalOutputFile(std::string path)
{
    InputScraper::SetTraversalOutputFile(std::move(path));
}

const std::string &
ParmParse::TraversalOutputFile()
{
    return InputScraper::TraversalOutputFile();
}

void
ParmParse::WriteInputTreeJson(std::ostream &os)
{
    InputScraper::WriteInputTreeJson(os);
}

void
ParmParse::WriteInputTreeJsonFile(const std::string &path)
{
    InputScraper::WriteInputTreeJsonFile(path);
}

bool
ParmParse::IgnoreInTraversalMode(   std::string note,
                                    const std::source_location &location)
{
    if (!InTraversalMode()) return false;
    InputScraper::RecordTraversalIgnore(location, std::move(note));
    return true;
}

void
ParmParse::PrintTraversalBranch(std::string name, const std::string &value)
{
    InputScraper::PrintTraversalBranch(*this, std::move(name), value);
}

void
ParmParse::RecordInput( std::string name,
                        std::string directive,
                        const std::source_location &location,
                        std::vector<std::string> options,
                        std::optional<std::string> default_value)
{
    InputScraper::RecordInput(  *this, std::move(name), std::move(directive),
                                location, std::move(options), std::move(default_value));
}

void
ParmParse::RecordConstraint(std::string kind,
                            int count,
                            std::vector<std::string> members,
                            std::vector<std::string> units,
                            const std::source_location &location)
{
    InputScraper::RecordConstraint( *this, std::move(kind), count, std::move(members),
                                    std::move(units), location);
}
}

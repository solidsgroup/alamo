#ifndef IO_INPUTSCRAPER_H
#define IO_INPUTSCRAPER_H

#include <optional>
#include <ostream>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "Set/Base.H"

namespace IO
{
class ParmParse;

class InputScraper
{
public:
    struct FlowNode;

    struct FlowBranch
    {
        std::string value;
        std::vector<FlowNode> children;
    };

    struct FlowNode
    {
        std::string kind = "input";
        std::string input;
        std::vector<FlowBranch> branches;
        std::vector<FlowNode> children;
    };

    struct InputNode
    {
        struct Condition
        {
            std::string path;
            std::string value;
        };

        struct Constraint
        {
            std::string kind;
            int count = 0;
            std::vector<std::string> members;
            std::vector<std::string> units;
            std::source_location location;
            std::vector<Condition> conditions;
        };

        std::string kind = "scope";
        std::string name;
        std::string full_name;
        std::string directive;
        std::source_location location;
        std::vector<std::string> options;
        bool required = false;
        bool has_default = false;
        bool has_default_value = false;
        std::string default_value;
        std::vector<std::vector<Condition>> contexts;
        std::vector<Constraint> constraints;
        std::vector<InputNode> children;
    };

    struct TraversalIgnore
    {
        std::string note;
        std::source_location location;
    };

    class TraversalBranchScope
    {
    public:
        TraversalBranchScope(ParmParse &pp, std::string name, std::string value);
        ~TraversalBranchScope();

    private:
        bool active = false;
    };

    static void SetTraversalMode(bool enabled);
    static bool InTraversalMode();
    static bool ShouldExecute();
    static const InputNode &InputTree();
    static void ClearInputTree();
    static void SetTraversalOutputFile(std::string path);
    static const std::string &TraversalOutputFile();
    static void WriteInputTreeJson(std::ostream &os);
    static void WriteInputTreeJsonFile(const std::string &path);

    static std::string DirectiveName(const std::source_location &location);

    static void RecordInput(ParmParse &pp,
                            std::string name,
                            std::string directive,
                            const std::source_location &location,
                            std::vector<std::string> options = {},
                            std::optional<std::string> default_value = std::nullopt);
    static void RecordConstraint(   ParmParse &pp,
                                    std::string kind,
                                    int count,
                                    std::vector<std::string> members,
                                    std::vector<std::string> units,
                                    const std::source_location &location);
    static void RecordTraversalIgnore(  const std::source_location &location,
                                        std::string note);
    static void CaptureSequenceTemplate(ParmParse &pp,
                                        const std::string &sequence_name,
                                        const std::string &template_name);
    static void PrintTraversalBranch(ParmParse &pp, std::string name, const std::string &value);

    template<typename T>
    static std::string InputValueString(const T &value)
    {
        std::ostringstream os;
        os << value;
        return os.str();
    }

    static std::string InputValueString(const std::string &value) { return value; }
    static std::string InputValueString(const char *value) { return std::string(value); }

    template<typename T>
    static std::string InputArrayString(const std::vector<T> &values)
    {
        std::ostringstream os;
        for (std::size_t i = 0; i < values.size(); i++)
        {
            if (i) os << " ";
            os << InputValueString(values[i]);
        }
        return os.str();
    }

    static std::string InputArrayString(const Set::Vector &value)
    {
        std::ostringstream os;
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            if (i) os << " ";
            os << value(i);
        }
        return os.str();
    }

    static std::string InputArrayString(const Set::Matrix &value)
    {
        std::ostringstream os;
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
            {
                if (i || j) os << " ";
                os << value(i, j);
            }
        }
        return os.str();
    }

private:
    static InputNode input_tree;
    static int traversal_print_depth;
    static bool traversal_mode;
    static std::string traversal_output_file;
    static std::vector<InputNode::Condition> traversal_conditions;
    static std::vector<FlowNode> input_flow;
    static std::vector<FlowNode> *traversal_flow;
    static std::vector<std::vector<FlowNode> *> traversal_flow_stack;
    static std::vector<TraversalIgnore> traversal_ignores;

    static InputNode &GetPath(InputNode &root, const std::string &path);
    static std::string KindForDirective(const std::string &directive,
                                        const std::vector<std::string> &options = {});
    static std::string MergeKind(const std::string &old_kind, const std::string &new_kind);
    static void AddContext(InputNode &node, const std::vector<InputNode::Condition> &conditions);
    static void RecordFlowInput(const std::string &path,
                                const std::string &directive,
                                const std::vector<std::string> &options);
    static void PushTraversalCondition(ParmParse &pp, std::string name, std::string value);
    static void PopTraversalCondition();
    static void PushTraversalBranch(ParmParse &pp, const std::string &name,
                                    const std::string &value);
    static void PopTraversalBranch();
};
}

#endif

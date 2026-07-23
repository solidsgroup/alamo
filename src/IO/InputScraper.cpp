#include "IO/InputScraper.H"
#include "IO/JSON.H"
#include "IO/ParmParse.H"

#include <fstream>
#include <iostream>
#include <utility>

#include "AMReX_ParallelDescriptor.H"

namespace
{
IO::InputScraper::InputNode &
GetChild(IO::InputScraper::InputNode &node, const std::string &name)
{
    for (auto &child : node.children)
        if (child.name == name)
            return child;

    IO::InputScraper::InputNode child;
    child.name = name;
    child.full_name = node.full_name.empty() ? name : node.full_name + "." + name;
    node.children.push_back(child);
    return node.children.back();
}

bool
SameConditions(const std::vector<IO::InputScraper::InputNode::Condition> &a,
               const std::vector<IO::InputScraper::InputNode::Condition> &b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
        if (a[i].path != b[i].path || a[i].value != b[i].value)
            return false;
    return true;
}

void
WriteConditions(std::ostream &os,
                const std::vector<IO::InputScraper::InputNode::Condition> &conditions,
                int indent)
{
    os << "[";
    for (std::size_t i = 0; i < conditions.size(); i++)
    {
        if (i) os << ",";
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        os << "{\"path\": \"" << IO::JSON::Escape(conditions[i].path)
           << "\", \"value\": \"" << IO::JSON::Escape(conditions[i].value) << "\"}";
    }
    if (!conditions.empty())
    {
        os << "\n";
        IO::JSON::Indent(os, indent);
    }
    os << "]";
}

void
WriteContextsField(std::ostream &os, bool &first, int indent,
                   const std::vector<std::vector<IO::InputScraper::InputNode::Condition>> &contexts)
{
    IO::JSON::Comma(os, first, indent);
    os << "\"contexts\": [";
    for (std::size_t i = 0; i < contexts.size(); i++)
    {
        if (i) os << ",";
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        WriteConditions(os, contexts[i], indent + 2);
    }
    if (!contexts.empty())
    {
        os << "\n";
        IO::JSON::Indent(os, indent);
    }
    os << "]";
}

void
WriteConstraint(std::ostream &os,
                const IO::InputScraper::InputNode::Constraint &constraint,
                int indent)
{
    bool first = true;
    os << "{";
    IO::JSON::WriteStringField(os, first, indent + 2, "kind", constraint.kind);
    IO::JSON::WriteIntField(os, first, indent + 2, "count", constraint.count);
    IO::JSON::WriteStringArrayField(os, first, indent + 2, "members", constraint.members);
    if (!constraint.units.empty())
        IO::JSON::WriteStringArrayField(os, first, indent + 2, "units", constraint.units);
    if (!constraint.conditions.empty())
    {
        IO::JSON::Comma(os, first, indent + 2);
        os << "\"conditions\": ";
        WriteConditions(os, constraint.conditions, indent + 2);
    }
    IO::JSON::WriteSourceField(os, first, indent + 2, constraint.location);
    os << "\n";
    IO::JSON::Indent(os, indent);
    os << "}";
}

void
WriteNode(std::ostream &os, const IO::InputScraper::InputNode &node, int indent)
{
    bool first = true;
    os << "{";
    IO::JSON::WriteStringField(os, first, indent + 2, "kind", node.kind);
    IO::JSON::WriteStringField(os, first, indent + 2, "name", node.name);
    IO::JSON::WriteStringField(os, first, indent + 2, "path", node.full_name);
    if (!node.directive.empty())
        IO::JSON::WriteStringField(os, first, indent + 2, "directive", node.directive);
    if (!node.options.empty())
        IO::JSON::WriteStringArrayField(os, first, indent + 2, "options", node.options);
    if (node.required)
        IO::JSON::WriteBoolField(os, first, indent + 2, "required", node.required);
    if (node.has_default)
    {
        IO::JSON::WriteBoolField(os, first, indent + 2, "has_default", node.has_default);
        if (node.has_default_value)
            IO::JSON::WriteStringField(os, first, indent + 2, "default_value", node.default_value);
    }
    if (!node.directive.empty())
        IO::JSON::WriteSourceField(os, first, indent + 2, node.location);
    if (!node.contexts.empty())
        WriteContextsField(os, first, indent + 2, node.contexts);
    if (!node.constraints.empty())
    {
        IO::JSON::Comma(os, first, indent + 2);
        os << "\"constraints\": [";
        for (std::size_t i = 0; i < node.constraints.size(); i++)
        {
            if (i) os << ",";
            os << "\n";
            IO::JSON::Indent(os, indent + 4);
            WriteConstraint(os, node.constraints[i], indent + 4);
        }
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        os << "]";
    }
    if (!node.children.empty())
    {
        IO::JSON::Comma(os, first, indent + 2);
        os << "\"children\": [";
        for (std::size_t i = 0; i < node.children.size(); i++)
        {
            if (i) os << ",";
            os << "\n";
            IO::JSON::Indent(os, indent + 4);
            WriteNode(os, node.children[i], indent + 4);
        }
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        os << "]";
    }
    os << "\n";
    IO::JSON::Indent(os, indent);
    os << "}";
}

void
WriteFlowNodes(std::ostream &os,
               const std::vector<IO::InputScraper::FlowNode> &nodes,
               int indent);

void
WriteFlowNode(std::ostream &os, const IO::InputScraper::FlowNode &node, int indent)
{
    bool first = true;
    os << "{";
    IO::JSON::WriteStringField(os, first, indent + 2, "kind", node.kind);
    IO::JSON::WriteStringField(os, first, indent + 2, "input", node.input);
    if (!node.branches.empty())
    {
        IO::JSON::Comma(os, first, indent + 2);
        os << "\"branches\": [";
        for (std::size_t i = 0; i < node.branches.size(); i++)
        {
            if (i) os << ",";
            os << "\n";
            IO::JSON::Indent(os, indent + 4);
            os << "{";
            bool branch_first = true;
            IO::JSON::WriteStringField(os, branch_first, indent + 6,
                                       "value", node.branches[i].value);
            IO::JSON::Comma(os, branch_first, indent + 6);
            os << "\"children\": ";
            WriteFlowNodes(os, node.branches[i].children, indent + 6);
            os << "\n";
            IO::JSON::Indent(os, indent + 4);
            os << "}";
        }
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        os << "]";
    }
    os << "\n";
    IO::JSON::Indent(os, indent);
    os << "}";
}

void
WriteFlowNodes(std::ostream &os,
               const std::vector<IO::InputScraper::FlowNode> &nodes,
               int indent)
{
    os << "[";
    for (std::size_t i = 0; i < nodes.size(); i++)
    {
        if (i) os << ",";
        os << "\n";
        IO::JSON::Indent(os, indent + 2);
        WriteFlowNode(os, nodes[i], indent + 2);
    }
    if (!nodes.empty())
    {
        os << "\n";
        IO::JSON::Indent(os, indent);
    }
    os << "]";
}
}

namespace IO
{
InputScraper::InputNode InputScraper::input_tree;
int InputScraper::traversal_print_depth = 0;
bool InputScraper::traversal_mode = false;
std::string InputScraper::traversal_output_file = "alamo-inputs.schema.json";
std::vector<InputScraper::InputNode::Condition> InputScraper::traversal_conditions;
std::vector<InputScraper::FlowNode> InputScraper::input_flow;
std::vector<InputScraper::FlowNode> *InputScraper::traversal_flow = &InputScraper::input_flow;
std::vector<std::vector<InputScraper::FlowNode> *> InputScraper::traversal_flow_stack;

InputScraper::InputNode &
InputScraper::GetPath(InputNode &root, const std::string &path)
{
    InputNode *node = &root;
    std::size_t start = 0;
    while (start < path.size())
    {
        std::size_t end = path.find('.', start);
        std::string part = path.substr(start, end - start);
        if (!part.empty())
            node = &GetChild(*node, part);
        if (end == std::string::npos)
            break;
        start = end + 1;
    }
    return *node;
}

void
InputScraper::AddContext(InputNode &node, const std::vector<InputNode::Condition> &conditions)
{
    for (const auto &context : node.contexts)
        if (SameConditions(context, conditions))
            return;
    node.contexts.push_back(conditions);
}

void
InputScraper::RecordFlowInput(const std::string &path,
                              const std::string &directive,
                              const std::vector<std::string> &options)
{
    if (!traversal_flow) traversal_flow = &input_flow;

    FlowNode *flow_node = nullptr;
    for (auto it = traversal_flow->rbegin(); it != traversal_flow->rend(); ++it)
    {
        if (it->input == path)
        {
            flow_node = &*it;
            break;
        }
    }

    if (!flow_node)
    {
        traversal_flow->push_back(FlowNode());
        flow_node = &traversal_flow->back();
        flow_node->input = path;
    }

    if (KindForDirective(directive) != "switch") return;

    flow_node->kind = "switch";
    for (const auto &option : options)
    {
        bool found = false;
        for (const auto &branch : flow_node->branches)
            if (branch.value == option)
                found = true;
        if (!found)
        {
            FlowBranch branch;
            branch.value = option;
            flow_node->branches.push_back(std::move(branch));
        }
    }
}

std::string
InputScraper::KindForDirective(const std::string &directive)
{
    if (directive == "query_switch" || directive == "query_if") return "switch";
    if (directive == "select" || directive == "select_default") return "switch";
    if (directive == "select_enumerate") return "sequence";
    if (directive == "query_enumerate" || directive == "queryarr_enumerate" || directive == "queryclass_enumerate") return "sequence";
    if (directive == "queryclass") return "scope";
    return "parameter";
}

std::string
InputScraper::MergeKind(const std::string &old_kind, const std::string &new_kind)
{
    if (old_kind == "scope") return new_kind;
    if (new_kind == "switch" || new_kind == "sequence") return new_kind;
    return old_kind;
}

std::string
InputScraper::DirectiveName(const std::source_location &location)
{
    std::string name = location.function_name();

    std::size_t paren = name.find('(');
    if (paren != std::string::npos)
        name = name.substr(0, paren);

    std::size_t scope = name.rfind("::");
    if (scope != std::string::npos)
        name = name.substr(scope + 2);

    std::size_t space = name.rfind(' ');
    if (space != std::string::npos)
        name = name.substr(space + 1);

    return name;
}

void
InputScraper::SetTraversalMode(bool enabled)
{
    traversal_mode = enabled;
    traversal_conditions.clear();
    if (enabled)
        ClearInputTree();
}

bool
InputScraper::InTraversalMode()
{
    return traversal_mode;
}

bool
InputScraper::ShouldExecute()
{
    return !traversal_mode;
}

const InputScraper::InputNode &
InputScraper::InputTree()
{
    return input_tree;
}

void
InputScraper::ClearInputTree()
{
    input_tree = InputNode();
    input_tree.kind = "scope";
    input_flow.clear();
    traversal_flow = &input_flow;
    traversal_flow_stack.clear();
}

void
InputScraper::SetTraversalOutputFile(std::string path)
{
    traversal_output_file = std::move(path);
}

const std::string &
InputScraper::TraversalOutputFile()
{
    return traversal_output_file;
}

void
InputScraper::WriteInputTreeJson(std::ostream &os)
{
    os << "{\n";
    os << "  \"schema_version\": 2,\n";
    os << "  \"format\": \"alamo.input_schema\",\n";
    os << "  \"root\": ";
    WriteNode(os, input_tree, 2);
    os << ",\n  \"flow\": ";
    WriteFlowNodes(os, input_flow, 2);
    os << "\n}\n";
}

void
InputScraper::WriteInputTreeJsonFile(const std::string &path)
{
    if (path.empty()) return;
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    std::ofstream os(path);
    WriteInputTreeJson(os);
}

void
InputScraper::PrintTraversalBranch(ParmParse &pp, std::string name, const std::string &value)
{
    if (!InTraversalMode()) return;
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    for (int i = 0; i < traversal_print_depth; i++)
        std::cout << "  ";

    std::cout << "-> " << pp.full(name) << " = " << value << std::endl;
}

void
InputScraper::RecordInput(ParmParse &pp,
                          std::string name,
                          std::string directive,
                          const std::source_location &location,
                          std::vector<std::string> options,
                          std::optional<std::string> default_value)
{
    if (!InTraversalMode()) return;

    InputNode &node = GetPath(input_tree, pp.full(name));
    node.kind = MergeKind(node.kind, KindForDirective(directive));
    node.directive = directive;
    node.location = location;
    if (!options.empty())
        node.options = std::move(options);
    node.required = node.required || directive.find("required") != std::string::npos;
    node.has_default = node.has_default || directive.find("default") != std::string::npos || default_value.has_value();
    if (default_value.has_value())
    {
        node.has_default_value = true;
        node.default_value = std::move(*default_value);
    }
    AddContext(node, traversal_conditions);
    RecordFlowInput(node.full_name, directive, node.options);

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    for (int i = 0; i < traversal_print_depth; i++)
        std::cout << "  ";

    std::cout << node.full_name << " [" << node.directive << "]";
    if (!node.options.empty())
    {
        std::cout << " {";
        for (std::size_t i = 0; i < node.options.size(); i++)
        {
            if (i) std::cout << ", ";
            std::cout << node.options[i];
        }
        std::cout << "}";
    }
    std::cout << std::endl;
}

void
InputScraper::RecordConstraint(ParmParse &pp,
                               std::string kind,
                               int count,
                               std::vector<std::string> members,
                               std::vector<std::string> units,
                               const std::source_location &location)
{
    if (!InTraversalMode()) return;

    for (auto &member : members)
        member = pp.full(member);

    InputNode &scope = GetPath(input_tree, pp.getPrefix());
    InputNode::Constraint constraint;
    constraint.kind = std::move(kind);
    constraint.count = count;
    constraint.members = std::move(members);
    constraint.units = std::move(units);
    constraint.location = location;
    constraint.conditions = traversal_conditions;
    scope.constraints.push_back(std::move(constraint));
}

void
InputScraper::PushTraversalCondition(ParmParse &pp, std::string name, std::string value)
{
    if (!InTraversalMode()) return;
    traversal_conditions.push_back({pp.full(name), std::move(value)});
}

void
InputScraper::PopTraversalCondition()
{
    if (!traversal_conditions.empty())
        traversal_conditions.pop_back();
}

void
InputScraper::PushTraversalBranch(ParmParse &pp, const std::string &name,
                                  const std::string &value)
{
    if (!traversal_flow) traversal_flow = &input_flow;

    const std::string path = pp.full(name);
    FlowNode *switch_node = nullptr;
    for (auto it = traversal_flow->rbegin(); it != traversal_flow->rend(); ++it)
    {
        if (it->input == path && it->kind == "switch")
        {
            switch_node = &*it;
            break;
        }
    }

    if (!switch_node)
    {
        traversal_flow->push_back(FlowNode());
        switch_node = &traversal_flow->back();
        switch_node->kind = "switch";
        switch_node->input = path;
    }

    FlowBranch *selected_branch = nullptr;
    for (auto &branch : switch_node->branches)
    {
        if (branch.value == value)
        {
            selected_branch = &branch;
            break;
        }
    }
    if (!selected_branch)
    {
        FlowBranch branch;
        branch.value = value;
        switch_node->branches.push_back(std::move(branch));
        selected_branch = &switch_node->branches.back();
    }

    traversal_flow_stack.push_back(traversal_flow);
    traversal_flow = &selected_branch->children;
}

void
InputScraper::PopTraversalBranch()
{
    if (traversal_flow_stack.empty())
    {
        traversal_flow = &input_flow;
        return;
    }
    traversal_flow = traversal_flow_stack.back();
    traversal_flow_stack.pop_back();
}

InputScraper::TraversalBranchScope::TraversalBranchScope(ParmParse &pp, std::string name, std::string value)
{
    if (InTraversalMode())
    {
        active = true;
        PushTraversalCondition(pp, name, value);
        PushTraversalBranch(pp, name, value);
        traversal_print_depth++;
    }
}

InputScraper::TraversalBranchScope::~TraversalBranchScope()
{
    if (active)
    {
        PopTraversalBranch();
        PopTraversalCondition();
        traversal_print_depth--;
    }
}
}

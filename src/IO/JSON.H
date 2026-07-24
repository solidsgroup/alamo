#ifndef IO_JSON_H
#define IO_JSON_H

#include <ostream>
#include <source_location>
#include <string>
#include <vector>

namespace IO
{
class JSON
{
public:
    static std::string Escape(const std::string &value)
    {
        std::string ret;
        for (char c : value)
        {
            switch (c)
            {
            case 92: ret += "\\\\"; break;
            case 34: ret += "\\\""; break;
            case 8: ret += "\\b"; break;
            case 12: ret += "\\f"; break;
            case 10: ret += "\\n"; break;
            case 13: ret += "\\r"; break;
            case 9: ret += "\\t"; break;
            default: ret += c; break;
            }
        }
        return ret;
    }

    static void Indent(std::ostream &os, int indent)
    {
        for (int i = 0; i < indent; i++) os << " ";
    }

    static void Comma(std::ostream &os, bool &first, int indent)
    {
        if (!first) os << ",\n";
        first = false;
        Indent(os, indent);
    }

    static void WriteStringField(   std::ostream &os, bool &first, int indent,
                                    const std::string &name, const std::string &value)
    {
        Comma(os, first, indent);
        os << "\"" << name << "\": \"" << Escape(value) << "\"";
    }

    static void WriteBoolField( std::ostream &os, bool &first, int indent,
                                const std::string &name, bool value)
    {
        Comma(os, first, indent);
        os << "\"" << name << "\": " << (value ? "true" : "false");
    }

    static void WriteIntField(  std::ostream &os, bool &first, int indent,
                                const std::string &name, int value)
    {
        Comma(os, first, indent);
        os << "\"" << name << "\": " << value;
    }

    static void WriteStringArray(std::ostream &os, const std::vector<std::string> &values)
    {
        os << "[";
        for (std::size_t i = 0; i < values.size(); i++)
        {
            if (i) os << ", ";
            os << "\"" << Escape(values[i]) << "\"";
        }
        os << "]";
    }

    static void WriteStringArrayField(  std::ostream &os, bool &first, int indent,
                                        const std::string &name,
                                        const std::vector<std::string> &values)
    {
        Comma(os, first, indent);
        os << "\"" << name << "\": ";
        WriteStringArray(os, values);
    }

    static void WriteSourceField(   std::ostream &os, bool &first, int indent,
                                    const std::source_location &location)
    {
        Comma(os, first, indent);
        os << "\"source\": {";
        os << "\"file\": \"" << Escape(location.file_name()) << "\", ";
        os << "\"line\": " << location.line() << ", ";
        os << "\"function\": \"" << Escape(location.function_name()) << "\"";
        os << "}";
    }
};
}

#endif

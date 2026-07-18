#ifndef MODEL_SPECIES_REGISTRY_H
#define MODEL_SPECIES_REGISTRY_H

#include <cctype>
#include <string>
#include <vector>

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Model
{
namespace Species
{
// Maps user-defined species identifiers to stable component indices and
// mechanics descriptors. Strings are resolved during parsing, not in kernels.
enum class Mechanics
{
    Fluid,
    DeformableSolid,
    RigidSolid
};

class Registry
{
public:
    void Parse(IO::ParmParse& pp)
    {
        pp.queryarr_required("species.names", names);
        if (names.empty())
            Util::Exception(INFO, "species.names must contain at least one identifier");

        mechanics.resize(names.size());
        for (int n = 0; n < Size(); ++n)
        {
            const std::string& id = names[n];
            if (!IsIdentifier(id))
                Util::Exception(INFO, id, " is not a valid species identifier");
            for (int m = 0; m < n; ++m)
                if (id == names[m])
                    Util::Exception(INFO, "Duplicate species identifier ", id);

            std::string type;
            pp.query_required(id + ".mechanics", type);
            if (type == "fluid")
                mechanics[n] = Mechanics::Fluid;
            else if (type == "deformable_solid")
                mechanics[n] = Mechanics::DeformableSolid;
            else if (type == "rigid_solid")
                mechanics[n] = Mechanics::RigidSolid;
            else
                Util::Exception(INFO, type, " is not a valid mechanics type for species ", id);
        }
    }

    static bool IsIdentifier(const std::string& id)
    {
        if (id.empty() ||
            !(std::isalpha(static_cast<unsigned char>(id[0])) || id[0] == '_'))
            return false;
        for (char c : id)
            if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_'))
                return false;
        return true;
    }

    int Size() const { return static_cast<int>(names.size()); }

    int Index(const std::string& id) const
    {
        for (int n = 0; n < Size(); ++n)
            if (names[n] == id) return n;
        return -1;
    }

    int Require(const std::string& id, const std::string& key) const
    {
        const int index = Index(id);
        if (index < 0)
            Util::Exception(INFO, "Unknown species identifier ", id, " in ", key);
        return index;
    }

    const std::string& Name(const int index) const { return names[index]; }
    const std::vector<std::string>& Names() const { return names; }
    Mechanics GetMechanics(const int index) const { return mechanics[index]; }

private:
    std::vector<std::string> names;
    std::vector<Mechanics> mechanics;
};
}
}

#endif

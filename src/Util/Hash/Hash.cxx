#include "Hash.H"

using namespace Util::Hash;

int main(int argc, char *argv[])
{
    std::filesystem::path relSrcPath("../../");
    std::filesystem::path exePath = std::filesystem::absolute(argv[0]).parent_path();
    std::filesystem::path srcPath = std::filesystem::canonical(exePath/relSrcPath);
    std::cout << srcPath << std::endl;
    getFinalHash(srcPath);
    return 0;
}

#include "Util/Zephyr.H"

#include <chrono>
#include <cstring>
#include <iostream>
#include <spawn.h>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

extern char** environ;

namespace
{
bool enabled = false;
pid_t sidecar_pid = -1;
}

namespace Util::Zephyr
{
void Enable()
{
    enabled = true;
}

void Start(const std::string& output_directory)
{
    if (!enabled || output_directory.empty() || sidecar_pid > 0) return;

    std::vector<std::string> arguments = {
        "zph",
        "watch",
        output_directory,
        "--pid",
        std::to_string(getpid()),
    };

    std::vector<char*> argv;
    argv.reserve(arguments.size() + 1);
    for (std::string& argument : arguments) argv.push_back(argument.data());
    argv.push_back(nullptr);

    const int error = posix_spawnp(
        &sidecar_pid,
        "zph",
        nullptr,
        nullptr,
        argv.data(),
        environ
    );
    if (error != 0)
    {
        sidecar_pid = -1;
        std::cerr << "Zephyr: unable to start zph: " << std::strerror(error)
                  << ". Install zph or run without --post." << std::endl;
        return;
    }

    std::cout << "Zephyr: posting " << output_directory << " with zph (PID "
              << sidecar_pid << ")" << std::endl;
}

void Stop()
{
    if (sidecar_pid <= 0) return;

    int status = 0;
    if (waitpid(sidecar_pid, &status, WNOHANG) == sidecar_pid)
    {
        sidecar_pid = -1;
        return;
    }

    // ALAMO has already written terminal metadata. Give zph a chance to
    // observe it and exit, but never signal the sidecar: for a very short run,
    // it may still be starting and a signal here can strand a remote run in
    // the "starting" state. If it needs longer, it will observe our PID exit.
    for (int attempt = 0; attempt < 30; ++attempt)
    {
        if (waitpid(sidecar_pid, &status, WNOHANG) == sidecar_pid)
        {
            sidecar_pid = -1;
            return;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}
}

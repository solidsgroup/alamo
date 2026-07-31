#include "IO/OutputLog.H"

#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <streambuf>
#include <string>

namespace
{

class LogSink
{
public:
    void Write(const char* data, std::streamsize size)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_state == State::Buffering)
            m_pending.append(data, static_cast<std::size_t>(size));
        else if (m_state == State::Open)
            m_file.write(data, size);
    }

    int Sync()
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_state != State::Open) return 0;
        m_file.flush();
        return m_file.good() ? 0 : -1;
    }

    void Open(const std::string& path, bool append)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_file.open(path, std::ios::out | (append ? std::ios::app : std::ios::trunc));
        if (!m_file.is_open())
            throw std::runtime_error("Could not open output log " + path);

        m_state = State::Open;
        m_file.write(
            m_pending.data(),
            static_cast<std::streamsize>(m_pending.size()));
        m_pending.clear();
        m_file.flush();
    }

    void Disable()
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_pending.clear();
        m_state = State::Disabled;
    }

    void Close()
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_file.is_open())
        {
            m_file.flush();
            m_file.close();
        }
        m_pending.clear();
        m_state = State::Disabled;
    }

private:
    enum class State
    {
        Buffering,
        Open,
        Disabled
    };

    std::mutex m_mutex;
    std::ofstream m_file;
    std::string m_pending;
    State m_state = State::Buffering;
};

class TeeBuffer
    : public std::streambuf
{
public:
    TeeBuffer(std::streambuf* terminal, LogSink& log)
        : m_terminal(terminal), m_log(log)
    {
    }

protected:
    int_type overflow(int_type character) override
    {
        if (traits_type::eq_int_type(character, traits_type::eof()))
            return traits_type::not_eof(character);

        const char value = traits_type::to_char_type(character);
        if (traits_type::eq_int_type(
                m_terminal->sputc(value), traits_type::eof()))
            return traits_type::eof();
        m_log.Write(&value, 1);
        return character;
    }

    std::streamsize xsputn(
        const char* data, std::streamsize size) override
    {
        const std::streamsize written = m_terminal->sputn(data, size);
        m_log.Write(data, written);
        return written;
    }

    int sync() override
    {
        const int terminal_status = m_terminal->pubsync();
        const int log_status = m_log.Sync();
        return terminal_status == 0 && log_status == 0 ? 0 : -1;
    }

private:
    std::streambuf* m_terminal;
    LogSink& m_log;
};

LogSink log_sink;
std::streambuf* terminal_stdout = std::cout.rdbuf();
std::streambuf* terminal_stderr = std::cerr.rdbuf();
std::streambuf* terminal_stdlog = std::clog.rdbuf();
TeeBuffer stdout_buffer(terminal_stdout, log_sink);
TeeBuffer stderr_buffer(terminal_stderr, log_sink);
TeeBuffer stdlog_buffer(terminal_stdlog, log_sink);
bool initialized = false;

}

namespace IO::OutputLog
{

void Initialize()
{
    if (initialized) return;
    std::cout.rdbuf(&stdout_buffer);
    std::cerr.rdbuf(&stderr_buffer);
    std::clog.rdbuf(&stdlog_buffer);
    initialized = true;
}

void Open(const std::string& path, bool append)
{
    log_sink.Open(path, append);
}

void DisableFile()
{
    log_sink.Disable();
}

void Finalize()
{
    if (!initialized) return;

    std::cout.flush();
    std::cerr.flush();
    std::clog.flush();
    std::cout.rdbuf(terminal_stdout);
    std::cerr.rdbuf(terminal_stderr);
    std::clog.rdbuf(terminal_stdlog);
    log_sink.Close();
    initialized = false;
}

}

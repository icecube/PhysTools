// pipe.cpp
// Work with named pipes with a C++ friendly interface


#include "PhysTools/pipe.h"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <csignal>
#include <fstream>

#include <fcntl.h>
#include <unistd.h>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <sys/wait.h>

#include <boost/version.hpp>

using namespace std;
namespace io = boost::iostreams;

namespace phys_tools {

const int       READ = 0;
const int       WRITE = 1;

Pipe::Pipe (ofstream& in, const string command)
{
    init (&in, command, NULL, NULL);
}

Pipe::Pipe (const string command, ifstream& out)
{
    init (NULL, command, &out, NULL);
}

Pipe::Pipe (const string command, ifstream& out, ifstream& err)
{
    init (NULL, command, &out, &err);
}

Pipe::Pipe (ofstream& in, const string command, ifstream& out)
{
    init (&in, command, &out, NULL);
}

Pipe::Pipe (ofstream& in, const string command, ifstream& out,
            ifstream& err)
{
    init (&in, command, &out, &err);
}

void
Pipe::init (ofstream* in,
            const string command,
            ifstream* out,
            ifstream* err)
{
    buf_in = 0;
    buf_out = 0;
    buf_err = 0;
    int     p_stdin[2];
    int     p_stdout[2];
    int     p_stderr[2];
    if (pipe (p_stdin) || pipe (p_stdout) || pipe (p_stderr)) {
        throw runtime_error ("could not open pipe");
    }
    pid = fork ();
    if (pid < 0) {
        throw runtime_error ("could not fork process");
    }
    if (pid == 0) {
        if (in) {
            close (p_stdin[WRITE]);
            dup2 (p_stdin[READ], READ);
        }
        if (out) {
            close (p_stdout[READ]);
            dup2 (p_stdout[WRITE], WRITE);
        }
        if (err) {
            close (p_stderr[READ]);
            dup2 (p_stderr[WRITE], 2);
        }
        execl ("/bin/sh", "sh", "-c", command.c_str (), NULL);
        throw runtime_error ("could not execute command");
    }
    else {
        assert (pid > 0);
        typedef io::stream_buffer<io::file_descriptor_sink> in_type;
        typedef io::stream_buffer<io::file_descriptor_source> out_type;
        if (in) {
			in->tie(NULL);
            close (p_stdin[READ]);
            fd_in = p_stdin[WRITE];
#if BOOST_VERSION >= 104400
            buf_in = new in_type (fd_in,io::never_close_handle);		
#else
            buf_in = new in_type (fd_in,false);
#endif
            ostream* in_ostream = in;
            in_ostream->rdbuf (buf_in);
        }
        if (out) {
			out->tie(in);
            close (p_stdout[WRITE]);
            fd_out = p_stdout[READ];
#if BOOST_VERSION >= 104400
            buf_out = new out_type (fd_out,io::never_close_handle);
#else
            buf_out = new out_type (fd_out,false);
#endif
            istream* out_istream = out;
            out_istream->rdbuf (buf_out);
        }
        if (err) {
			err->tie(in);
            close (p_stderr[WRITE]);
            fd_err = p_stderr[READ];
#if BOOST_VERSION >= 104400
            buf_err = new out_type (fd_err,io::never_close_handle);
#else
            buf_err = new out_type (fd_err,false);
#endif
            istream* err_istream = err;
            err_istream->rdbuf (buf_err);
        }
    }
}

Pipe::~Pipe ()
{
    // it looks like if we don't do this, we leak resources which,
    // and evenutally we're not allowed to open any more pipes!
    finish_input ();
    finish_output ();
}


void
Pipe::finish_input ()
{
    if (buf_in) {
        close (fd_in);
    }
}


void
Pipe::finish_output ()
{
    if (buf_out) {
        close (fd_out);
    }
    if (buf_err) {
        close (fd_err);
    }
}


void
Pipe::terminate ()
{
    kill (pid, SIGTERM);
}


void
Pipe::wait ()
{
    waitpid (pid, NULL, 0);
}

} // namespace phys_tools

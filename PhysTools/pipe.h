/// @file pipe.hpp
/// @brief Work with named pipes with a C++ friendly interface
///
/// The provided implementation of the Pipe class declared in this
/// file requires the Boost.IOStreams library to files known by file
/// descriptors to standard c++ fstreams.


#ifndef PIPE_HPP
#define PIPE_HPP 1

#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>


namespace phys_tools {

/// @defgroup System System Interface
/// @brief A simple, C++-friendly interface to other binaries.
///
/// This set of classes provides an interface to access some
/// functionality provided by other programs (e.g. Gnuplot).  The root
/// of its functionality is the Pipe class.  In principle, Pipe can be
/// used to leverage any external binary that provides a well-defined
/// command line interface.
///
/// @{

/// @brief A named pipe class.
///
/// This Pipe class tries to make it very easy to setup pipes to and
/// from a command.  The constructors make it intuitive to start a
/// command, optionally sending it information on its stdin via a
/// std::ofstream, and optionally collecting information from it on its
/// stdout or stderr via std::ifstreams.
///
/// When the Pipe is destroyed, the file descriptors corresponding to
/// the application's ends of the Pipe are closed (otherwise, there is
/// a "file descriptor leak").  Therefore, the lifetime of the Pipe
/// determines the useful lifetime of the associated std::fstreams.
class Pipe {
public:

    /// @brief Pipe data to a command.
    /// @param in       An ofstream.
    /// @param command  The command.
    ///
    /// This constructor sets up a one-directional pipe to the given
    /// command.
    Pipe (std::ofstream& in, const std::string command);

    /// @brief Pipe data from a command.
    /// @param command  The command.
    /// @param out      An ifstream.
    ///
    /// This constructor sets up a one-directional pipe from the given
    /// command.
    Pipe (const std::string command, std::ifstream& out);

    /// @brief Pipe data from a command.
    /// @param command  The command.
    /// @param out      An ifstream.
    /// @param err      An ifstream.
    ///
    /// This constructor sets up a pipe from the given command, where
    /// its stdout output is directed into @c out and its stderr output
    /// is directed into @c err.
    Pipe (const std::string command, std::ifstream& out,
          std::ifstream& err);
    
    /// @brief Pipe data to and from a command.
    /// @param in       An ofstream.
    /// @param command  The command.
    /// @param out      An ifstream.
    ///
    /// This constructor sets up a bidirectional pipe to and from the
    /// given command, where data can be sent to it via @c in and read
    /// from its stdout via @c out.
    Pipe (std::ofstream& in, const std::string command, std::ifstream& out);

    /// @brief Pipe data to and from a command.
    /// @param in       An ofstream.
    /// @param command  The command.
    /// @param out      An ifstream.
    /// @param err      An ifstream.
    ///
    /// This constructor sets up a bidirectional pipe to and from the
    /// given command, where data can be sent to it via @c in, read
    /// from its stdout via @c out, and read from its stderr via @c err.
    Pipe (std::ofstream& in, const std::string command, std::ifstream& out,
          std::ifstream& err);

    ~Pipe ();

    /// @brief If this pipe has in input connection, close it.
    void            finish_input ();

    /// @brief If this pipe has an output connection (or connections)
    /// close it (or them)
    void            finish_output ();

    /// @brief Terminate the forked process.
    void            terminate ();

    /// @brief Wait for the forked process to terminate.
    void            wait ();

private:
    void            init (std::ofstream* in,
                          const std::string command,
                          std::ifstream* out,
                          std::ifstream* err);
    pid_t           pid;
    int             fd_in;
    int             fd_out;
    int             fd_err;
    std::streambuf* buf_in;
    std::streambuf* buf_out;
    std::streambuf* buf_err;
};

/// @}

} // namespace phys_tools

#endif  /* PIPE_HPP */

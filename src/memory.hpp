#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <cstdio>
#include <mutex>
#include <thread>

#include <mpi.h>

class MPI_guard
{
private:
    int size;
    int rank;
    std::string hostname;
public:
    MPI_guard(int *argc, char **argv[]);

    ~MPI_guard();
};


#endif // MEMORY_HPP

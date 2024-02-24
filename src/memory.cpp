#include <memory.hpp>

static std::once_flag init;
static std::once_flag dropped;

MPI_guard::MPI_guard(int *argc, char **argv[])
{
    std::call_once(init, [argc, argv](){
        MPI_Init(argc, argv);

        int _size, _rank;
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        printf("<MPI_guard:%d/%d> Init once! \n", _rank, _size);
    });

    // Get the world size (number of processes)
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the processor name (may not be available on all systems)
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    hostname.assign(processor_name, name_len);
}

MPI_guard::~MPI_guard()
{
    std::call_once(dropped, []() {
        int _size, _rank;
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        printf("<MPI_guard:%d/%d> Dropped once! \n", _rank, _size);
        MPI_Finalize();
    });
}
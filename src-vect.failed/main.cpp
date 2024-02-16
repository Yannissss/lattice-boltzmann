/*
** Code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the Bhatnagar-Gross-Krook collision step.
**
** The 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** A 2D grid:
**
**           cols
**       --- --- ---
**      | D | E | F |
** rows  --- --- ---
**      | A | B | C |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1D array:
**
**  --- --- --- --- --- ---
** | A | B | C | D | E | F |
**  --- --- --- --- --- ---
**
** Grid indicies are:
**
**          ny
**          ^       cols(ii)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(jj) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** Note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   ./d2q9-bgk input.params obstacles.dat
**
** Be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>

#include <driver.hpp>
#include <kernel.hpp>

#define FINALSTATEFILE "final_state.dat"
#define INITIALSTATEFILE "initial_state.dat"
#define AVVELSFILE "av_vels.dat"

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle
 * densities */
int initialise(const char *paramfile, const char *obstaclefile, t_param *params,
               t_speed *cells_ptr, t_speed *tmp_cells_ptr, int **obstacles_ptr,
               float **av_vels_ptr);

int write_values(const t_param params, t_speed cells, int *obstacles,
                 float *av_vels);

/* finalise, including freeing up allocated memory */
int finalise(const t_param *params, t_speed *cells_ptr, t_speed *tmp_cells_ptr,
             int **obstacles_ptr, float **av_vels_ptr);

/* utility functions */
void die(const char *message, const int line, const char *file);
void usage(const char *exe);

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char *argv[])
{
    char *paramfile = NULL;    /* name of the input parameter file */
    char *obstaclefile = NULL; /* name of a the input obstacle file */
    t_param params;            /* struct to hold parameter values */
    t_speed cells;             /* grid containing fluid densities */
    t_speed tmp_cells;         /* scratch space */
    int *obstacles = NULL;     /* grid indicating which cells are blocked */
    float *av_vels =
        NULL; /* a record of the av. velocity computed for each timestep */
    struct timeval timstr; /* structure to hold elapsed time */
    double tot_tic, tot_toc, init_tic, init_toc, comp_tic, comp_toc, col_tic,
        col_toc; /* floating point numbers to calculate elapsed wallclock time
                  */

    /* parse the command line */
    if (argc != 3)
    {
        usage(argv[0]);
    }
    else
    {
        paramfile = argv[1];
        obstaclefile = argv[2];
    }

    /* Total/init time starts here: initialise our data structures and load
     * values from file */
    gettimeofday(&timstr, NULL);
    tot_tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
    init_tic = tot_tic;
    initialise(paramfile, obstaclefile, &params, &cells, &tmp_cells, &obstacles,
               &av_vels);
#ifdef DEBUG
    write_values(params, cells, obstacles, av_vels);
#endif

    /* Init time stops here, compute time starts*/
    gettimeofday(&timstr, NULL);
    init_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
    comp_tic = init_toc;

    // BEGIN : Compute section
    lattice_boltzmann(params, cells, tmp_cells, obstacles, av_vels);
    // END : Compute section

    /* Compute time stops here, collate time starts*/
    gettimeofday(&timstr, NULL);
    comp_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
    col_tic = comp_toc;

    // Collate data from ranks here

    /* Total/collate time stops here.*/
    gettimeofday(&timstr, NULL);
    col_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
    tot_toc = col_toc;

    /* write final values and free memory */
    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n",
           calc_reynolds(params, cells, obstacles));
    printf("Elapsed Init time:\t\t\t%.6lf (s)\n", init_toc - init_tic);
    printf("Elapsed Compute time:\t\t\t%.6lf (s)\n", comp_toc - comp_tic);
    printf("Elapsed Collate time:\t\t\t%.6lf (s)\n", col_toc - col_tic);
    printf("Elapsed Total time:\t\t\t%.6lf (s)\n", tot_toc - tot_tic);
    write_values(params, cells, obstacles, av_vels);
    finalise(&params, &cells, &tmp_cells, &obstacles, &av_vels);

    return EXIT_SUCCESS;
}

int initialise(const char *paramfile, const char *obstaclefile, t_param *params,
               t_speed *cells_ptr, t_speed *tmp_cells_ptr, int **obstacles_ptr,
               float **av_vels_ptr)
{
    char message[1024]; /* message buffer */
    FILE *fp;           /* file pointer */
    int xx, yy;         /* generic array indices */
    int blocked;        /* indicates whether a cell is blocked by an obstacle */
    int retval;         /* to hold return value for checking */

    /* open the parameter file */
    fp = fopen(paramfile, "r");

    if (fp == NULL)
    {
        sprintf(message, "could not open input parameter file: %s", paramfile);
        die(message, __LINE__, __FILE__);
    }

    /* read in the parameter values */
    retval = fscanf(fp, "%d\n", &(params->nx));

    if (retval != 1)
        die("could not read param file: nx", __LINE__, __FILE__);

    retval = fscanf(fp, "%d\n", &(params->ny));

    if (retval != 1)
        die("could not read param file: ny", __LINE__, __FILE__);

    retval = fscanf(fp, "%d\n", &(params->maxIters));

    if (retval != 1)
        die("could not read param file: maxIters", __LINE__, __FILE__);

    retval = fscanf(fp, "%d\n", &(params->reynolds_dim));

    if (retval != 1)
        die("could not read param file: reynolds_dim", __LINE__, __FILE__);

    retval = fscanf(fp, "%f\n", &(params->density));

    if (retval != 1)
        die("could not read param file: density", __LINE__, __FILE__);

    retval = fscanf(fp, "%f\n", &(params->accel));

    if (retval != 1)
        die("could not read param file: accel", __LINE__, __FILE__);

    retval = fscanf(fp, "%f\n", &(params->omega));

    if (retval != 1)
        die("could not read param file: omega", __LINE__, __FILE__);

    /* and close up the file */
    fclose(fp);

    /*
    ** Allocate memory.
    **
    ** Remember C is pass-by-value, so we need to
    ** pass pointers into the initialise function.
    **
    ** NB we are allocating a 1D array, so that the
    ** memory will be contiguous.  We still want to
    ** index this memory as if it were a (row major
    ** ordered) 2D array, however.  We will perform
    ** some arithmetic using the row and column
    ** coordinates, inside the square brackets, when
    ** we want to access elements of this array.
    **
    ** Note also that we are using a structure to
    ** hold an array of 'speeds'.  We will allocate
    ** a 1D array of these structs.
    */

    float *ptr = NULL;

    /* main grid */
    for (int kk = 0; kk < NSPEEDS; kk++)
    {
        ptr = (float *)malloc(sizeof(float) * (params->ny * params->nx));

        if (ptr == NULL)
            die("cannot allocate memory for cells", __LINE__, __FILE__);

        cells_ptr->speeds[kk] = ptr;
    }

    /* 'helper' grid, used as scratch space */
    for (int kk = 0; kk < NSPEEDS; kk++)
    {
        ptr = (float *)malloc(sizeof(float) * (params->ny * params->nx));

        if (ptr == NULL)
            die("cannot allocate memory for tmp_cells", __LINE__, __FILE__);

        tmp_cells_ptr->speeds[kk] = ptr;
    }

    /* the map of obstacles */
    *obstacles_ptr = (int *)malloc(sizeof(int) * (params->ny * params->nx));

    if (*obstacles_ptr == NULL)
        die("cannot allocate column memory for obstacles", __LINE__, __FILE__);

    /* initialise densities */
    float w0 = params->density * 4.f / 9.f;
    float w1 = params->density / 9.f;
    float w2 = params->density / 36.f;

    for (int jj = 0; jj < params->ny; jj++)
    {
        for (int ii = 0; ii < params->nx; ii++)
        {
            /* centre */
            cells_ptr->speeds[0][ii + jj * params->nx] = w0;
            /* axis directions */
            cells_ptr->speeds[1][ii + jj * params->nx] = w1;
            cells_ptr->speeds[2][ii + jj * params->nx] = w1;
            cells_ptr->speeds[3][ii + jj * params->nx] = w1;
            cells_ptr->speeds[4][ii + jj * params->nx] = w1;
            /* diagonals */
            cells_ptr->speeds[5][ii + jj * params->nx] = w2;
            cells_ptr->speeds[6][ii + jj * params->nx] = w2;
            cells_ptr->speeds[7][ii + jj * params->nx] = w2;
            cells_ptr->speeds[8][ii + jj * params->nx] = w2;
        }
    }

    /* first set all cells in obstacle array to zero */
    for (int jj = 0; jj < params->ny; jj++)
    {
        for (int ii = 0; ii < params->nx; ii++)
        {
            (*obstacles_ptr)[ii + jj * params->nx] = 0;
        }
    }

    /* open the obstacle data file */
    fp = fopen(obstaclefile, "r");

    if (fp == NULL)
    {
        sprintf(message, "could not open input obstacles file: %s",
                obstaclefile);
        die(message, __LINE__, __FILE__);
    }

    /* read-in the blocked cells list */
    while ((retval = fscanf(fp, "%d %d %d\n", &xx, &yy, &blocked)) != EOF)
    {
        /* some checks */
        if (retval != 3)
            die("expected 3 values per line in obstacle file", __LINE__,
                __FILE__);

        if (xx < 0 || xx > params->nx - 1)
            die("obstacle x-coord out of range", __LINE__, __FILE__);

        if (yy < 0 || yy > params->ny - 1)
            die("obstacle y-coord out of range", __LINE__, __FILE__);

        if (blocked != 1)
            die("obstacle blocked value should be 1", __LINE__, __FILE__);

        /* assign to array */
        (*obstacles_ptr)[xx + yy * params->nx] = blocked;
    }

    /* and close the file */
    fclose(fp);

    /*
    ** allocate space to hold a record of the avarage velocities computed
    ** at each timestep
    */
    *av_vels_ptr = (float *)malloc(sizeof(float) * params->maxIters);

    return EXIT_SUCCESS;
}

int finalise(const t_param *params, t_speed *cells_ptr, t_speed *tmp_cells_ptr,
             int **obstacles_ptr, float **av_vels_ptr)
{
    /*
    ** free up allocated memory
    */
    for (int kk = 0; kk < NSPEEDS; kk++)
    {
        free(cells_ptr->speeds[kk]);
        cells_ptr->speeds[kk] = NULL;

        free(tmp_cells_ptr->speeds[kk]);
        tmp_cells_ptr->speeds[kk] = NULL;
    }

    free(*obstacles_ptr);
    *obstacles_ptr = NULL;

    free(*av_vels_ptr);
    *av_vels_ptr = NULL;

    return EXIT_SUCCESS;
}

int write_values(const t_param params, t_speed cells, int *obstacles,
                 float *av_vels)
{
    FILE *fp;                     /* file pointer */
    const float c_sq = 1.f / 3.f; /* sq. of speed of sound */
    float local_density;          /* per grid cell sum of densities */
    float pressure;               /* fluid pressure in grid cell */
    float u_x;                    /* x-component of velocity in grid cell */
    float u_y;                    /* y-component of velocity in grid cell */
    float u; /* norm--root of summed squares--of u_x and u_y */

#ifdef DEBUG
    /* Hack degueux pour créer une image initiale si build en mode debug */
    static int is_final = 0;
    if (!is_final)
    {
        fp = fopen(INITIALSTATEFILE, "w");
        is_final = !is_final;
    }
    else
#endif
        fp = fopen(FINALSTATEFILE, "w");

    if (fp == NULL)
    {
        die("could not open file output file", __LINE__, __FILE__);
    }

    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* an occupied cell */
            if (obstacles[ii + jj * params.nx])
            {
                u_x = u_y = u = 0.f;
                pressure = params.density * c_sq;
            }
            /* no obstacle */
            else
            {
                local_density = 0.f;

                for (int kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells.speeds[kk][ii + jj * params.nx];
                }

                /* compute x velocity component */
                u_x = (cells.speeds[1][ii + jj * params.nx] +
                       cells.speeds[5][ii + jj * params.nx] +
                       cells.speeds[8][ii + jj * params.nx] -
                       (cells.speeds[3][ii + jj * params.nx] +
                        cells.speeds[6][ii + jj * params.nx] +
                        cells.speeds[7][ii + jj * params.nx])) /
                      local_density;
                /* compute y velocity component */
                u_y = (cells.speeds[2][ii + jj * params.nx] +
                       cells.speeds[5][ii + jj * params.nx] +
                       cells.speeds[6][ii + jj * params.nx] -
                       (cells.speeds[4][ii + jj * params.nx] +
                        cells.speeds[7][ii + jj * params.nx] +
                        cells.speeds[8][ii + jj * params.nx])) /
                      local_density;
                /* compute norm of velocity */
                u = sqrtf((u_x * u_x) + (u_y * u_y));
                /* compute pressure */
                pressure = local_density * c_sq;
            }

            /* write to file */
            fprintf(fp, "%d %d %.12E %.12E %.12E %.12E %d\n", ii, jj, u_x, u_y,
                    u, pressure, obstacles[ii + params.nx * jj]);
        }
    }

    fclose(fp);

    fp = fopen(AVVELSFILE, "w");

    if (fp == NULL)
    {
        die("could not open file output file", __LINE__, __FILE__);
    }

    for (int ii = 0; ii < params.maxIters; ii++)
    {
        fprintf(fp, "%d:\t%.12E\n", ii, av_vels[ii]);
    }

    fclose(fp);

    return EXIT_SUCCESS;
}

void die(const char *message, const int line, const char *file)
{
    fprintf(stderr, "Error at line %d of file %s:\n", line, file);
    fprintf(stderr, "%s\n", message);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void usage(const char *exe)
{
    fprintf(stderr, "Usage: %s <paramfile> <obstaclefile>\n", exe);
    exit(EXIT_FAILURE);
}
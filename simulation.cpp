#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <time.h>
using namespace std;
using namespace Eigen;

double Simulation::run()
{
    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i= 0; i < temporal.NT; i++)
    {
        iteration(i);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    return elapsed;
}

void Simulation::iteration(int i)
{
    //periodic for now
    
    //grid.save_field_values
    /* grid.apply_bc(i); */
    for (Species species: list_species)
    {
        species.periodic_interpolate_fields(grid);
        species.velocity_push();
        grid.gather_charge(species);
        grid.gather_current(species);
    }
    grid.solve();
    for (Species species: list_species)
    {
        species.position_push();
        /* species.save_particle_values(i); */
        species.periodic_apply_bc(grid);
    }
}


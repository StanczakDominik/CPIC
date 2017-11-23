#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <time.h>
using namespace std;
using namespace Eigen;

Simulation::Simulation(Temporal& _temporal, Grid& _grid, string _filename, Species &species)
    : temporal(_temporal), grid(_grid), filename(_filename)
{
    list_species = std::vector<Species>();
    list_species.push_back(species);
}

Simulation::Simulation(Temporal& _temporal, Grid &_grid, string _filename, Species& species, Species& species2)
    : temporal(_temporal), grid(_grid), filename(_filename)
{
    list_species = std::vector<Species>();
    list_species.push_back(species);
    list_species.push_back(species2);
}

double Simulation::run()
{
    struct timespec start, finish;
    double elapsed;

    grid.initial_solve((bool)0);

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
    //grid.save_field_values
    /* grid.apply_bc(i); */
    for (Species species: list_species)
    {
        species.interpolate_fields(grid);
        species.velocity_push();
        species.gather_charge(grid);
        species.gather_current(grid);
    }
    grid.solve();
    for (Species species: list_species)
    {
        species.position_push();
        /* species.save_particle_values(i); */
        grid.apply_particle_bc(species);
    }
}


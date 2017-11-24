#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <time.h>
#include "H5Cpp.h"
using namespace std;
using namespace Eigen;
using namespace H5;

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
    grid.apply_bc(temporal.dt*i);
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
        grid.apply_particle_bc(species);
    }
}

void Simulation::save()
{
    H5File file(filename, H5F_ACC_TRUNC);
    hsize_t scalar_grid[1];
    scalar_grid[0] = grid.NG;

    hsize_t transversal_grid[2];
    transversal_grid[0] = grid.NG;
    transversal_grid[1] = 2;

    hsize_t vector_grid[2];
    vector_grid[0] = grid.NG;
    vector_grid[1] = 3;

    DataSpace scalar_dataspace(1, scalar_grid);
    DataSpace transversal_dataspace(2, transversal_grid);
    DataSpace vector_dataspace(2, vector_grid);

    DataSet grid_x = file.createDataSet(H5std_string("grid_x"), PredType::NATIVE_DOUBLE, scalar_dataspace);
    grid_x.write(grid.x.data(), PredType::NATIVE_DOUBLE);

    DataSet charge_density = file.createDataSet(H5std_string("grid_charge_density"), PredType::NATIVE_DOUBLE, scalar_dataspace);
    charge_density.write(grid.charge_density.data(), PredType::NATIVE_DOUBLE);

    DataSet current_density_x = file.createDataSet(H5std_string("grid_current_density_x"), PredType::NATIVE_DOUBLE, scalar_dataspace);
    current_density_x.write(grid.current_density_x.data(), PredType::NATIVE_DOUBLE);

    DataSet current_density_yz = file.createDataSet(H5std_string("grid_current_density_yz"), PredType::NATIVE_DOUBLE, transversal_dataspace);
    current_density_yz.write(grid.current_density_yz.data(), PredType::NATIVE_DOUBLE);

    DataSet electric_field = file.createDataSet(H5std_string("grid_electric_field"), PredType::NATIVE_DOUBLE, vector_dataspace);
    electric_field.write(grid.electric_field.data(), PredType::NATIVE_DOUBLE);

    DataSet magnetic_field = file.createDataSet(H5std_string("grid_magnetic_field"), PredType::NATIVE_DOUBLE, vector_dataspace);
    magnetic_field.write(grid.magnetic_field.data(), PredType::NATIVE_DOUBLE);
}

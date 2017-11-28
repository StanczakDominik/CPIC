#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <time.h>
using namespace std;
using namespace Eigen;

Simulation::Simulation(Temporal& _temporal, Grid& _grid, string _filename, Species* species)
    : temporal(_temporal), grid(_grid), filename(_filename)
{
    list_species = std::vector<Species *>();
    list_species.push_back(species);
}

Simulation::Simulation(Temporal& _temporal, Grid &_grid, string _filename, Species* species, Species* species2)
    : temporal(_temporal), grid(_grid), filename(_filename)
{
    list_species = std::vector<Species *>();
    list_species.push_back(species);
    list_species.push_back(species2);
}

void Simulation::grid_species_init()
{
    grid.apply_bc(0);
    grid.current_density_x = ArrayXd::Zero(grid.NG+3);
    grid.current_density_yz = ArrayX2d::Zero(grid.NG+4, 2);
    grid.charge_density = ArrayXd::Zero(grid.NG+1);
    for (int i = 0; i<(int)list_species.size(); i++)
    {
        list_species[i]->interpolate_fields(grid);
        list_species[i]->velocity_push(grid);
        list_species[i]->gather_charge(grid);
        list_species[i]->gather_current(grid);
        /* cout << "position push" << endl; */
        list_species[i]->position_push();
        /* cout << "apply particle bc" << endl; */
        list_species[i]->apply_particle_bc(grid);
    }
}


double Simulation::run()
{
    struct timespec start, finish;
    double elapsed;

    grid_species_init();

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i= 0; i < temporal.NT; i++)
    {
        /* printf("\rIteration%5d/%5d", i, temporal.NT); */
        iteration(i);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    return elapsed;
}

void Simulation::iteration(int i)
{
    cout << "\riteration " << i << ", apply bc";
    grid.apply_bc(temporal.dt*i);
    grid.current_density_x = ArrayXd::Zero(grid.NG+3);
    grid.current_density_yz = ArrayX2d::Zero(grid.NG+4, 2);
    grid.charge_density = ArrayXd::Zero(grid.NG+1);
    for (int j = 0; j<(int)list_species.size(); j++)
    {
        cout << "\riteration " << i << ", interpolate fields" ;
        list_species[j]->interpolate_fields(grid);
        cout << "\riteration " << i << ", velocity push" ;
        list_species[j]->velocity_push(grid);
        cout << "\riteration " << i << ", gather charge" ;
        list_species[j]->gather_charge(grid);
        cout << "\riteration " << i << ", gather current" ;
        list_species[j]->gather_current(grid);
        cout << "\riteration " << i << ", position push" ;
        list_species[j]->position_push();
        cout << "\riteration " << i << ", apply particle bc" ;
        list_species[j]->apply_particle_bc(grid);
    }
    cout << "\riteration " << i << ", grid solve" ;
    grid.solve();
}

void Simulation::save()
{

    /* H5File file(filename, H5F_ACC_TRUNC); */
    /* hsize_t scalar_grid[1]; */
    /* scalar_grid[0] = grid.NG; */

    /* hsize_t transversal_grid[2]; */
    /* transversal_grid[0] = grid.NG; */
    /* transversal_grid[1] = 2; */

    /* hsize_t vector_grid[2]; */
    /* vector_grid[0] = grid.NG; */
    /* vector_grid[1] = 3; */

    /* Group group = Group( file.createGroup( "/grid" )); */

    /* DataSpace scalar_dataspace(1, scalar_grid); */
    /* DataSpace transversal_dataspace(2, transversal_grid); */
    /* DataSpace vector_dataspace(2, vector_grid); */

    /* DataSet grid_x = file.createDataSet(H5std_string("x"), PredType::NATIVE_DOUBLE, scalar_dataspace); */
    /* grid_x.write(grid.x.data(), PredType::NATIVE_DOUBLE); */

    /* DataSet charge_density = group.createDataSet(H5std_string("rho"), PredType::NATIVE_DOUBLE, scalar_dataspace); */
    /* charge_density.write(grid.charge_density.data(), PredType::NATIVE_DOUBLE); */

    /* DataSet current_density_x = group.createDataSet(H5std_string("current_x"), PredType::NATIVE_DOUBLE, scalar_dataspace); */
    /* current_density_x.write(grid.current_density_x.data(), PredType::NATIVE_DOUBLE); */

    /* DataSet current_density_yz = group.createDataSet(H5std_string("current_yz"), PredType::NATIVE_DOUBLE, transversal_dataspace); */
    /* current_density_yz.write(grid.current_density_yz.data(), PredType::NATIVE_DOUBLE); */

    /* DataSet electric_field = group.createDataSet(H5std_string("Efield"), PredType::NATIVE_DOUBLE, vector_dataspace); */
    /* electric_field.write(grid.electric_field.data(), PredType::NATIVE_DOUBLE); */

    /* DataSet magnetic_field = group.createDataSet(H5std_string("Bfield"), PredType::NATIVE_DOUBLE, vector_dataspace); */
    /* magnetic_field.write(grid.magnetic_field.data(), PredType::NATIVE_DOUBLE); */
}

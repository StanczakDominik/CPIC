#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include "simulation.hpp"
using namespace std;
using namespace Eigen;

void Simulation::iteration(int i)
{
    //periodic for now
    
    //grid.save_field_values
    grid.apply_bc(i);
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
        species.apply_bc();
    }
}


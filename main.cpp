#include <iostream>
#include "grid.hpp"
#include "species.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <vector>
#include "H5Cpp.h"

float laser_wavelength = 1.064e-6;
float laser_intensity = 1e21;
float impulse_duration = 1e-13;

float length = 1.0655e-5 ;
float total_time = 2e-13 ;
float spatial_step = 7.7325e-9 ;

float moat_length_left_side = 3.093e-6 ;

double scaling = 9.847700361687114e+24;

string category_name = "benchmark_run";

float epsilon_zero = 8.854187817e-12;
float electric_charge = 1.60217662e-19;
float lightspeed = 299792458;
float proton_mass = 1.672621898e-27;
float electron_rest_mass = 9.10938356e-31;

double test_run(int n_macroparticles, int n_cells)
{
   cout << "Initializing variables" << endl;
   Temporal temp(spatial_step/lightspeed, total_time); 
   NonPeriodicSpecies electrons(n_macroparticles, -electric_charge, electron_rest_mass, scaling, temp.dt);
   NonPeriodicSpecies protons(n_macroparticles, electric_charge, proton_mass, scaling, temp.dt);
   NonPeriodicGrid g(n_cells, length, lightspeed, epsilon_zero, temp, laser_wavelength, laser_intensity, total_time/2.0, impulse_duration, 6);
   std::ostringstream stringStream;
   stringStream << n_macroparticles << "_" << n_cells << ".hdf5";
   std::string filename = stringStream.str();
   Simulation sim(temp, g, filename, &electrons, &protons);
   for (int i = 0; i<(int)sim.list_species.size(); i++)
   {
      sim.list_species[i]->distribute_uniformly(g, 0, moat_length_left_side, moat_length_left_side);
      sim.list_species[i]->apply_particle_bc(g);
   }

   cout << "Running sim" << endl;
   double runtime = sim.run();
   cout << "Running sim took " << runtime << " seconds" << endl;
   sim.save();
   return runtime;
}

int main()
{
   int n_particles[] = {100, 1000, 10000, 50000, 75000};
   int n_grid[] = {100, 500, 1000, 2000};
   for(int j = 0; j < 5; j++)
   {
     for (int i = 0; i < 4; i++)
     {
        int number_grid = n_grid[i];
        int number_particles = n_particles[j];
        if (number_particles >= number_grid)
           cout << number_particles << "," << number_grid << "," << test_run(number_particles, number_grid) << endl;
     }
   }
}


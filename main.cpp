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

float laser_wavelength = 1.064e-6
float laser_intensity = 1e21
float impulse_duration = 1e-13

float length = 1.0655e-5 # meters
float total_time = 2e-13 # seconds
float spatial_step = 7.7325e-9 # meters

float moat_length_left_side = 3.093e-6 # meters

double scaling = 9.847700361687114e+24

category_name = "benchmark_run"

epsilon_zero = 8.854187817e-12
electric_charge = 1.60217662e-19
lightspeed = 299792458
proton_mass = 1.672621898e-27
electron_rest_mass = 9.10938356e-31

double test_run(int n_macroparticles, int n_cells)
{
   cout << "Initializing variables" << endl;
   Temporal temp(spatial_step/lightspeed, total_time); 
   NonPeriodicSpecies electrons(n_macroparticles, -electric_charge, electron_rest_mass, scaling, temp.dt);
   NonPeriodicSpecies protons(n_macroparticles, electric_charge, proton_mass, scaling, temp.dt);
   NonPeriodicGrid g(n_cells, length, lightspeed, epsilon_0, temp, laser_wavelength, laser_intensity, total_time/2.0, impulse_duration, 6);
   Simulation sim(temp, g, "filename.hdf5", &e, &p);
   for (int i = 0; i<(int)sim.list_species.size(); i++)
   {
      sim.list_species[i]->distribute_uniformly(g, 0, moat_length_left_side, moat_length_left_side);
   }

   cout << "Running sim" << endl;
   double runtime = sim.run();
   cout << "Running sim took " << runtime << " seconds" << endl;
   cout << s.N_alive << endl;
   sim.save();
   return runtime;
}

int main()
{
   cout << test_run(10000, 1000) << endl;
}


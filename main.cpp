#include <iostream>
#include "grid.hpp"
#include "species.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <Eigen/Dense>
#include <fstream>
using namespace std;
using namespace Eigen;
#include <vector>

double laser_wavelength = 1.064e-6;
double laser_intensity = 1e21;
double impulse_duration = 1e-13;

double length = 1.0655e-5 ;
double total_time = 2e-13 ;
double spatial_step = 7.7325e-9 ;

double moat_length_left_side = 3.093e-6 ;

double scaling = 9.847700361687114e+24;

string category_name = "benchmark_run";

double epsilon_zero = 8.854187817e-12;
double electric_charge = 1.60217662e-19;
double lightspeed = 299792458;
double proton_mass = 1.672621898e-27;
double electron_rest_mass = 9.10938356e-31;

double test_run(int n_macroparticles, int n_cells, bool save)
{
   Temporal temp(n_cells, total_time); 
   cout << "Running for " << n_cells << " cells, " << n_macroparticles << " particles, " << temp.NT << " iterations " << endl;
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
      sim.list_species[i]->v = ArrayX3d::Constant(n_macroparticles,3, (1-2*i)*lightspeed/10000.0);
   }

   double runtime = sim.run(save);
   cout << "\rRunning sim took " << runtime << " seconds" << endl;
   return runtime;
}

int main()
{
   int n_particles[] = {100, 200, 500, 750, 1000, 1750, 2000, 2500, 5000, 10000, 20000, 50000};
   vector<int> nparticles(n_particles, n_particles + sizeof(n_particles) / sizeof(n_particles[0]));
   /* int n_particles[] = {1000}; */
   int number_grid = 1000;
   std::ofstream out("dane.csv");
   for(int number_particles: nparticles)
   {
     double runtime = test_run(number_particles, number_grid, false);
     out << number_particles << "," << runtime << endl;
   }
   out.close();
}


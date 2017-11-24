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


int main()
{
   cout << "Initializing variables" << endl;
   Temporal temp(2, 1.0); 
   NonPeriodicSpecies s(100, 1, 1, 1, temp.dt);
   s.v = ArrayX3d::Constant(s.N_alive, 3, 0.1);
   /* NonPeriodicGrid g(32, 1, 1, 1, temp, 1, 1, 0.5, 0.25, 2); */
   Grid g(32, 1, 1, 1, temp);
   s.distribute_uniformly(g, 1e-10, 0, 0);
   s.sinusoidal_position_perturbation(1e-3, 1, g);
   g.apply_particle_bc(s);
   Simulation sim(temp, g, "filename.hdf5", s);
   cout << "Running sim" << endl;
   double runtime = sim.run();
   cout << "Running sim took " << runtime << " seconds" << endl;
   cout << s.N_alive << endl;
   sim.save();
}

#include <iostream>
#include "grid.hpp"
#include "species.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <vector>

int main()
{
   cout << "Initializing variables" << endl;
   Temporal temp(1000, 1.0); 
   Species s(1000, 1, 1, 1);
   Grid g(32, 1, 1, 1, temp);
   s.distribute_uniformly(g, 1e-10, 0, 0);
   s.sinusoidal_position_perturbation(1e-3, 1, g);
   string file = "filename";
   Simulation sim(temp, g, file, s);
   cout << "Running sim" << endl;
   double runtime = sim.run();
   cout << "Running sim took " << runtime << " seconds" << endl;
   cout << sim.grid.charge_density << endl;
   cout << sim.list_species.at(0).x << endl;
}

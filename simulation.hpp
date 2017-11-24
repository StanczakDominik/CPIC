#ifndef SIMULATION_H
#define SIMULATION_H
#include <Eigen/Dense>
#include <vector>
#include "grid.hpp"
#include "species.hpp"
#include "temporal.hpp"
#include "H5Cpp.h"
using namespace std;
using namespace Eigen;

class Simulation
{
    public:
        Temporal temporal;

        Grid grid;
        std::vector<Species> list_species;
        
        string filename;
        /* string title; */
    Simulation(Temporal &temporal, Grid &grid, string filename, Species &species);
    Simulation(Temporal &temporal, Grid &grid, string filename, Species &species, Species &species2);

    void iteration(int i);
    double run();
};

#endif /* SIMULATION_H */

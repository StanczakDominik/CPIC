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
using namespace H5;

class Simulation
{
    public:
        Temporal temporal;

        Grid grid;
        std::vector<Species *> list_species;
        
        H5std_string filename;
        /* string title; */
    Simulation(Temporal &temporal, Grid &grid, string filename, Species *species);
    Simulation(Temporal &temporal, Grid &grid, string filename, Species *species, Species *species2);

    void iteration(int i);
    double run(bool save);
    void grid_species_init();
    void save(int i);
};

#endif /* SIMULATION_H */

#ifndef SIMULATION_H
#define SIMULATION_H
#include <Eigen/Dense>
#include <vector>
#include "temporal.hpp"
using namespace std;
using namespace Eigen;

class Grid;
class Species;
class Simulation
{
    public:
        Temporal temporal;

        Grid grid;
        std::vector<Species> list_species;
        
        string filename;
        string title;

    void iteration(int i);
    double run();
};

#endif /* SIMULATION_H */

#ifndef SIMULATION_H
#define SIMULATION_H
#include <Eigen/Dense>
#include "temporal.hpp"
#include <vector>
using namespace std;
using namespace Eigen;

class Grid;
class Species;
class Simulation
{
    public:
        Temporal timing;

        Grid grid;
        std::vector<Species> list_species;
        
        string filename;
        string title;

    void iteration(int i);
    void run();
};

#endif /* SIMULATION_H */

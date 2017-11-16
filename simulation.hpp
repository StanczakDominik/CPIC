#ifndef SIMULATION_H
#define SIMULATION_H
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Grid;
class Species;
class Simulation
{
    public:
        int NT;
        float T;
        float dt;
        ArrayXd t;

        Grid grid;
        // TODO list_species
        Species species;
        
        str filename; // TODO get str
        str title;

    void iteration(int i);
    void run();
};

#endif /* SIMULATION_H */

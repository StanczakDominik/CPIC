#include <Eigen/Dense>
#include "species.hpp"
#ifndef GRID_H
#define GRID_H
using namespace std;
using namespace Eigen;

class Species;
class Grid
{
    public:
        float c;
        float epsilon_0;
        int NG;
        float L;

        ArrayXd x;
        float dx;
        ArrayXd charge_density;
        ArrayXd current_density_x;
        ArrayX2d current_density_yz;
        ArrayX3d electric_field;
        ArrayX3d magnetic_field;
    Grid(int _NG, float _L, float _c, float _epsilon_0);
    void initial_solve(bool neutralize);
    void solve();
    ArrayXd bincount(ArrayXd cell_numbers, ArrayXd weights, int minlength);
    void gather_charge(Species s);
    void gather_charge_periodic(Species s);
    void gather_current(Species s);
    void gather_current_periodic(Species s);
};

#endif /* GRID_H */

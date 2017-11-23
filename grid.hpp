#ifndef GRID_H
#define GRID_H
#include <Eigen/Dense>
#include "temporal.hpp"
using namespace std;
using namespace Eigen;

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
        Temporal temporal;
    Grid(int _NG, float _L, float _c, float _epsilon_0, Temporal &_temporal);
    void initial_solve(bool neutralize);
    void solve();
    ArrayXd bincount(const Ref<ArrayXd>& cell_numbers, const Ref<ArrayXd>& weights, int minlength);
    void apply_bc(int iteration);
};

class NonPeriodicGrid : public Grid
{
    void apply_bc(int iteration);
};

#endif /* GRID_H */


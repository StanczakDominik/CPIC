#ifndef GRID_H
#define GRID_H
#include <Eigen/Dense>
#include "temporal.hpp"

using namespace std;
using namespace Eigen;

class Grid
{
    public:
        double c;
        double epsilon_0;
        int NG;
        double L;

        ArrayXd x;
        double dx;
        ArrayXd charge_density;
        ArrayXd current_density_x;
        ArrayX2d current_density_yz;
        ArrayX3d electric_field;
        ArrayX3d magnetic_field;
        Temporal temporal;
    Grid(int _NG, double _L, double _c, double _epsilon_0, Temporal &_temporal);
    virtual ~Grid();
    void initial_solve(bool neutralize);
    void solve();
    ArrayXd bincount(const Ref<ArrayXd>& cell_numbers, const Ref<ArrayXd>& weights, int minlength);
    void apply_bc(double t);
};

class NonPeriodicGrid : public Grid
{
    public:
        double laser_omega;
        double laser_amplitude;
        double envelope_center_t;
        double envelope_width;
        double envelope_power;
        NonPeriodicGrid(int _NG, double _L, double _c, double _epsilon_0,
                Temporal &_temporal, double _laser_wavelength, double
                _laser_intensity, double _envelope_center_t, double
                _envelope_width, double _envelope_power);
        void apply_bc(double t);
    private:
        double _taui;
        double _tau;
        double _t_0;
};

#endif /* GRID_H */


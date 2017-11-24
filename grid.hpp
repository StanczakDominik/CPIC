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
    void apply_bc(float t);
};

class NonPeriodicGrid : public Grid
{
    public:
        float laser_omega;
        float laser_amplitude;
        float envelope_center_t;
        float envelope_width;
        float envelope_power;
        NonPeriodicGrid(int _NG, float _L, float _c, float _epsilon_0, Temporal &_temporal,
            float _laser_wavelength, float _laser_intensity, float _envelope_center_t, float _envelope_width,
            float _envelope_power);
        void apply_bc(float t);
    private:
        float _taui;
        float _tau;
        float _t_0;
};

#endif /* GRID_H */


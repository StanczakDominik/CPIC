#include <iostream>
#include <Eigen/Dense>
#include "grid.hpp"
#include "species.hpp"
using namespace std;
using namespace Eigen;

ArrayXd PythonLinSpaced(int NG, float x_min, float x_max)
{
    // Equivalent to python's numpy.linspace(... endpoint = False;
    ArrayXd x(NG);
    float dx = (x_max-x_min)/NG;
    for (int i = 1; i < NG; i++)
    {
        x(i) = i * dx;
    }
    return x;
}


Grid::Grid(int _NG, float _L, float _c, float _epsilon_0)
{
    // compute effective charges and masses of macroparticles
    L = _L;
    NG = _NG;
    c = _c;
    epsilon_0 = _epsilon_0;
    // allocate position and velocity arrays
    x = PythonLinSpaced(NG, 0, L);
    dx = x(1) - x(0);
}

ArrayXd Grid::gather_density(Species s)
{
    ArrayXd logical_coordinates = floor(s.x / dx);
    ArrayXd charge_to_right = (s.x / dx) - logical_coordinates;
    ArrayXd charge_to_left = 1.0 - charge_to_right;
    /* charge_hist_to_right = np.bincount(logical_coordinates+1, charge_to_right, minlength=x.size+1) */
    /* charge_hist_to_left = np.bincount(logical_coordinates, charge_to_left, minlength=x.size+1) */
    /* return charge_hist_to_right + charge_hist_to_left */
    ArrayXd charge_density(NG);
    return charge_density;
}

void Grid::initial_solve(bool neutralize)
{
    /* rho_F = fft(rho); */
    /* if(neutralize) */
    /* { */
    /*     rho_F(0) = 0; */
    /* } */
    /* field_F = rho_F / (1j * k * epsilon) */
    /* return fft.ifft(field_F).real) */
}

void Grid::solve()
{
    ArrayXd Fplus = 0.5 * (electric_field.col(0) + c * magnetic_field.col(1));
    ArrayXd Fminus = 0.5 * (electric_field.col(0) - c * magnetic_field.col(1));
    ArrayXd Gplus = 0.5 * (electric_field.col(1) + c * magnetic_field.col(0));
    ArrayXd Gminus = 0.5 * (electric_field.col(1) - c * magnetic_field.col(0));

    /* Fplus.tail(NG-1) = Fplus.head(NG-1) - 0.5 * dt * current.block(2, NG) / epsilon_0; */
    /* Gplus.tail(NG-1) = Gplus.head(NG-1) - 0.5 * dt * current.block(2, NG) / epsilon_0; */
    /* Fminus.head(NG-1) = Fminus.tail(NG-1) - 0.5 * dt * current.block(2, NG) / epsilon_0; */
    /* Gminus.head(NG-1) = Gminus.tail(NG-1) - 0.5 * dt * current.block(2, NG) / epsilon_0; */
}


void test_grid()
{
    Grid g = Grid(10, 1, 1, 1);

    cout << g.x << endl;
    cout << g.x.mean() << endl;
    cout << g.c << endl;
    cout << g.epsilon_0 << endl;
    cout << g.x.size() << endl;
}


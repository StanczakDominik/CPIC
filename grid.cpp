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
    ArrayXd charge_density(NG);
    return charge_density;
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


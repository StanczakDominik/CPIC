#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
using namespace std;
using namespace Eigen;

Species::Species(int _N, float _q, float _m, float _scaling)
{
    // compute effective charges and masses of macroparticles
    N = _N;
    q = _q;
    m = _m;
    scaling = _scaling;
    eff_q = q * scaling;
    eff_m = m * scaling;
    // allocate position and velocity arrays
    x = ArrayXd(N);
    v = ArrayX3d(N, 3);
}

void Species::position_push()
{
    x += v.col(0) * dt;
}





void test_species()
{
    Species s = Species(3, 1, 1, 10.0);
    s.x << 1.5, 3.5, 6.5;

    cout << s.x.mean() << endl;
    cout << s.q << endl;
    cout << s.v.size() << endl;
}


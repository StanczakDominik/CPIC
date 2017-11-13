#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
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


void Species::velocity_push()
{
    float c = 1; // todo reference from grid
    v.colwise() /= (1 - v.pow(2).rowwise().sum()/pow(c,2)).sqrt();
    ArrayX3d half_force = (q * 0.5 / m * dt) * E;
    v += half_force;

    ArrayX3d t = B * q * dt / (2 * m * (1 + v.pow(2).rowwise().sum()/pow(c,2)).sqrt());
    cout << v.rows() << endl << v.cols() << endl;
    cout << t.rows() << endl << t.cols() << endl;
    /* ArrayX3d uprime = v + v.rowwise().cross(t.rowwise()); */
    /* ArrayXd divisor = 1 - (v/c).pow(2).rowwise().sum(); */
    /* cout << divisor.rows() << ", " << divisor.cols() << endl; */
    /* cout << v.rows() << ", " << v.cols() << endl; */
    /* v.colwise() /= divisor; */
    
    return;
}

void interpolate_fields(Grid g)
{
    /* ArrayXd logical_coordinates = floor(s.x / dx); */
    /* ArrayXd right_fractions  = (s.x / dx) - logical_coordinates; */
    
    /* ArrayXd charge_hist_to_right = bincount(logical_coordinates+1, charge_to_right, NG+1); */
    /* ArrayXd charge_hist_to_left = bincount(logical_coordinates, 1-charge_to_right, NG+1); */
    /* ArrayXd charge_density(NG); */
    /* charge_density = charge_hist_to_right + charge_hist_to_left; */
    // TODO
    return ;
}



void test_species()
{
    Species s = Species(3, 1, 1, 10.0);
    s.x << 1.5, 3.5, 6.5;

    cout << s.x.mean() << endl;
    cout << s.q << endl;
    cout << s.v.size() << endl;
}


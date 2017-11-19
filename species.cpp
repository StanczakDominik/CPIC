#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include <cmath>
using namespace std;
using namespace Eigen;

Species::Species(int _N, float _q, float _m, float _scaling)
{
    // compute effective charges and masses of macroparticles
    N = _N;
    q = _q;
    m = _m;
    N_alive = N;
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


double Species::velocity_push()
{
    float c = 1; // todo reference from grid
    v.colwise() /= (1 - v.pow(2).rowwise().sum()/pow(c,2)).sqrt();
    ArrayX3d half_force = (q * 0.5 / m * dt) * E;
    v += half_force;

    ArrayX3d t = B * q * dt / (2 * m);
    t.colwise() /= (1 + v.pow(2).rowwise().sum()/pow(c,2)).sqrt();
    ArrayX3d t2 = 2 * t;
    t2.colwise() /= 1 + t.pow(2).rowwise().sum();

    MatrixX3d uprime = v.matrix();
    for (int i = 0; i < N_alive; i++)
    {
        uprime.row(i) += v.row(i).matrix().cross(t.row(i).matrix());
        v.row(i).matrix() += uprime.row(i).matrix().cross(t2.row(i).matrix());
    }
    v += half_force;

    ArrayXd final_gamma =(1 + v.pow(2).rowwise().sum()/pow(c,2)).sqrt();
    v.colwise() /= final_gamma;

    return (final_gamma - 1).sum() * eff_m * c * c;
}

void Species::periodic_interpolate_fields(Grid g)
{
    ArrayXd logical_coordinates = floor(x / g.dx);
    ArrayXd right_fractions  = (x / g.dx) - logical_coordinates;

    for (int i =0; i < N_alive; i ++)
    {
        E(i) = (1 - right_fractions(i)) * g.electric_field(i+1) + right_fractions(i) * g.electric_field((i+1) % g.NG + 1);
        B(i) = (1 - right_fractions(i)) * g.magnetic_field(i+1) + right_fractions(i) * g.magnetic_field((i+1) % g.NG + 1);
    }
}


void Species::aperiodic_interpolate_fields(Grid g)
{
    ArrayXd logical_coordinates = floor(x / g.dx);
    ArrayXd right_fractions  = (x / g.dx) - logical_coordinates;

    for (int i =0; i < N_alive; i ++)
    {
        E(i) = (1 - right_fractions(i)) * g.electric_field(i+1) + right_fractions(i) * E(i+2);
        B(i) = (1 - right_fractions(i)) * g.magnetic_field(i+1) + right_fractions(i) * B(i+2);
    }
}

void Species::periodic_apply_bc(Grid g)
{
    x = x - (g.L * (x/g.L).floor());
}

void Species::aperiodic_apply_bc(Grid g)
{
    ArrayXi indices = ((0 < x) * (x < g.L)).cast<int>();
    int N_alive_new = indices.sum();
    if (N_alive != N_alive_new)
    {
        ArrayXd new_x(N_alive_new);
        ArrayX3d new_v(N_alive_new, 3);
        int j = 0;
        for (int i = 0;  i < N_alive; i++)
        {
            if (indices(i))
            {
                new_x(j) = x(i);
                new_v.row(j) = v.row(i);
                j++;
            }
        }
        N_alive = N_alive_new;
        x = new_x;
    }

    /* alive = (0 <= x) & (x < g.L) */
    /* if N_alive: */
    /*     x = x[alive] */
    /*     v = v[alive] */
    /* self.N_alive = alive.sum() */
}

void test_species()
{
    Species s = Species(3, 1, 1, 10.0);
    s.x << 1.5, 3.5, 6.5;

    cout << s.x.mean() << endl;
    cout << s.q << endl;
    cout << s.v.size() << endl;
}

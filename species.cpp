#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include <cmath>
using namespace std;
using namespace Eigen;

Species::Species(int _N, float _q, float _m, float _scaling)
    : q(_q), m(_m), N(_N), N_alive(_N), scaling(_scaling), eff_q(q*scaling),
    eff_m(m*scaling), x(N), v(N, 3), E(N, 3), B(N, 3)
{
    // compute effective charges and masses of macroparticles
    /* N = _N; */
    /* q = _q; */
    /* m = _m; */
    /* N_alive = N; */
    /* scaling = _scaling; */
    /* eff_q = q * scaling; */
    /* eff_m = m * scaling; */
    /* // allocate position and velocity arrays */
    /* x = ArrayXd(N); */
    /* v = ArrayX3d(N, 3); */
}

void Species::position_push()
{
    x += v.col(0) * dt;
}

void Species::distribute_uniformly(Grid &g, float shift, float start_moat, float end_moat)
{
    x = ArrayXd::LinSpaced(N,
            start_moat + g.L / N * 1e-10,
            g.L-end_moat); // TODO check endpoint, python has false
    x += shift * N / g.L / 10.0;
    apply_bc(g);
}

void Species::sinusoidal_position_perturbation(float amplitude, int mode, Grid &g)
{
    x += amplitude * cos(2*mode*M_PI*x / g.L);
    apply_bc(g);
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

void Species::interpolate_fields(Grid &g)
{
    for (int i =0; i < N_alive; i ++)
    {
        int on_grid = (x(i) / g.dx);
        double right_fraction = (x(i) / g.dx) - on_grid;
        E.row(i) = (1 - right_fraction) * g.electric_field.row(on_grid+1) + right_fraction * g.electric_field.row((on_grid+1) % g.NG + 1);
        B.row(i) = (1 - right_fraction) * g.magnetic_field.row(on_grid+1) + right_fraction * g.magnetic_field.row((on_grid+1) % g.NG + 1);
    }
}


void NonPeriodicSpecies::interpolate_fields(Grid &g)
{
    for (int i =0; i < N_alive; i ++)
    {
        int on_grid = (x(i) / g.dx);
        double right_fraction = (x(i) / g.dx) - on_grid;
        E.row(i) = (1 - right_fraction) * g.electric_field.row(on_grid+1) + right_fraction * g.electric_field.row(on_grid+2);
        B.row(i) = (1 - right_fraction) * g.magnetic_field.row(on_grid+1) + right_fraction * g.magnetic_field.row(on_grid+2);
    }
}

void Species::apply_bc(Grid &g)
{
    x = x - (g.L * (x/g.L).floor());
}

void NonPeriodicSpecies::apply_bc(Grid &g)
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
}

// CHARGE AND CURRENT DEPOSITION

void Species::gather_charge_computation(Grid &g)
{
    ArrayXd logical_coordinates = floor(x / g.dx);
    ArrayXd charge_to_right = 1.0 - ((x / g.dx) - logical_coordinates);
   
    ArrayXd charge_hist_to_left = g.bincount(logical_coordinates, charge_to_right, g.NG+1);
    logical_coordinates += 1;
    charge_to_right = 1 - charge_to_right;
    ArrayXd charge_hist_to_right = g.bincount(logical_coordinates, charge_to_right, g.NG+1);

    g.charge_density = charge_hist_to_right + charge_hist_to_left;
}

void Species::gather_charge(Grid &g)
{
    gather_charge_computation(g);
    g.charge_density.head(1) += g.charge_density.tail(1);
}

void NonPeriodicSpecies::gather_charge(Grid &g)
{
    gather_charge_computation(g);
}

typedef Array<bool,Dynamic,1> ArrayXb;

void Species::gather_current_computation(Grid &g)
{
    float epsilon = g.dx * 1e-10;
    for (int i=0; i < N_alive; i++)
    {
        float x_velocity = v(i,0);
        bool active = v.row(i).any();
        float time_left = dt;
        float xp = x(i);
        while(active)
        {
            int logical_coordinate = (int)floor(xp/g.dx);

            bool particle_in_left_half = x(i) / g.dx - logical_coordinate <= 0.5;
            float s, t1;
            if (particle_in_left_half)
            {
                if(x_velocity > 0)
                {
                    t1 = -(xp - logical_coordinate * g.dx) / x_velocity;
                    s = logical_coordinate * g.dx - epsilon;
                }
                else
                {
                    t1 = ((logical_coordinate + 0.5) * g.dx - xp) / x_velocity;
                    s = (logical_coordinate + 0.5) * g.dx + epsilon;
                }
            }
            else // particle in right half
            {
                if(x_velocity > 0)
                {
                    t1 = ((logical_coordinate + 1 ) * g.dx)/ x_velocity;
                    s = (logical_coordinate + 1) * g.dx + epsilon;
                }
                else
                {
                    t1 = -(xp - (logical_coordinate + 0.5) * g.dx) / x_velocity;
                    s = (logical_coordinate + 0.5) * g.dx + epsilon;
                }
            }

            float time_overflow = time_left - t1;
            bool switches_cells = time_overflow  > 0;
            float time_in_this_iteration = switches_cells ? t1 : time_left;
            time_in_this_iteration = (x_velocity == 0) ? dt : time_in_this_iteration;

            int logical_coordinate_long = particle_in_left_half ? logical_coordinate: logical_coordinate + 1;
            int logical_coordinate_trans = particle_in_left_half ? logical_coordinate-1: logical_coordinate + 1;

            int sign = (int)(particle_in_left_half) * 2 - 1;
            float distance_to_center = (logical_coordinate + 0.5) * g.dx - xp;
            float s0 = 1 - sign * distance_to_center / g.dx;
            float change_in_coverage = sign * x_velocity * time_in_this_iteration / g.dx;
            float s1 = s0 + change_in_coverage;
            float w = 0.5 * (s0 + s1);

            Array3d j_contribution = v.row(i) * eff_q / dt * time_in_this_iteration;
            float y_contribution_to_current_cell = w * j_contribution(1);
            float z_contribution_to_current_cell = w * j_contribution(2);
            float y_contribution_to_next_cell = (1-w) * j_contribution(1);
            float z_contribution_to_next_cell = (1-w) * j_contribution(2);

            g.current_density_x(logical_coordinate_long + 1) += j_contribution(0);
            g.current_density_yz(logical_coordinate + 2, 0) += y_contribution_to_current_cell;
            g.current_density_yz(logical_coordinate + 2, 1) += z_contribution_to_current_cell;
            g.current_density_yz(logical_coordinate_trans + 2, 0) += y_contribution_to_next_cell;
            g.current_density_yz(logical_coordinate_trans + 2, 1) += z_contribution_to_next_cell;

            active = switches_cells;
            time_left = time_overflow;
            xp = s;
        }
    }
}

void Species::gather_current(Grid &g)
{
    gather_current_computation(g);
    g.current_density_yz.block(g.current_density_x.rows() - 4,0 , 2, 2) += g.current_density_yz.topRows(2);
    g.current_density_yz.block(2, 0 , 2, 2) += g.current_density_yz.bottomRows(2);
    g.current_density_x(g.current_density_x.rows() -3) += g.current_density_x(0); // TODO check from end
    g.current_density_x(0) = 0;
    g.current_density_x.segment(2, 1) += g.current_density_x.tail(2);
    g.current_density_x.tail(2) = 0;
}

void NonPeriodicSpecies::gather_current(Grid &g)
{
    gather_current_computation(g);
    g.current_density_yz.topRows(2) = 0;
    g.current_density_yz.bottomRows(2) = 0;
    g.current_density_x(0) = 0;
    g.current_density_x.tail(2) = 0;

}


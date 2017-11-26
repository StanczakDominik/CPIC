#include <iostream>
#include <Eigen/Dense>
#include "species.hpp"
#include "grid.hpp"
#include <cmath>
using namespace std;
using namespace Eigen;

Species::Species(int _N, float _q, float _m, float _scaling, float _dt)
    : q(_q), m(_m), N(_N), N_alive(_N), scaling(_scaling), eff_q(q*scaling),
    eff_m(m*scaling), dt(_dt), E(N, 3), B(N, 3)
{
    x = ArrayXd::Zero(N);
    v = ArrayX3d::Zero(N,3);
    E = ArrayX3d::Zero(N,3);
    B = ArrayX3d::Zero(N,3);
}

Species::~Species()
{
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
    apply_particle_bc(g);
}

void Species::sinusoidal_position_perturbation(float amplitude, int mode, Grid &g)
{
    x += amplitude * cos(2*mode*M_PI*x / g.L);
    apply_particle_bc(g);
}



double Species::velocity_push(Grid& g)
{
    double c = g.c;
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


// CHARGE AND CURRENT DEPOSITION

void Species::gather_charge_computation(Grid &g)
{
    ArrayXd logical_coordinates = floor(x / g.dx);
    ArrayXd charge_to_right = (x / g.dx) - logical_coordinates;
   
    g.charge_density = ArrayXd::Zero(g.NG+1);
    for (int i = 0; i < N_alive; i++)
    {
        int logical_coordinate = floor(x(i) / g.dx);
        double charge_to_right = x(i) / g.dx - logical_coordinate;
        g.charge_density[logical_coordinate] += 1-charge_to_right;
        g.charge_density[logical_coordinate+1] += charge_to_right;
    }
}

void Species::gather_charge(Grid &g)
{
    gather_charge_computation(g);
    g.charge_density.head(1) += g.charge_density.tail(1);
}

typedef Array<bool,Dynamic,1> ArrayXb;

void Species::gather_current_computation(Grid &g)
{
    double eps = g.dx * 1e-4;
    for (int i=0; i < N_alive; i++)
    {
        double x_velocity = v(i,0);
        bool active = v.row(i).any();
        double time_left = dt;
        double xp = x(i);
        /* int emergency_counter = 0; */
        while(active)
        {
            /* emergency_counter++; */
            int logical_coordinate = (int)floor(xp/g.dx);

            bool particle_in_left_half = xp / g.dx - logical_coordinate <= 0.5;
            double s, t1;
            if (particle_in_left_half)
            {
                if(x_velocity < 0) // case 1
                {
                    t1 = -(xp - logical_coordinate * g.dx) / x_velocity;
                    s = logical_coordinate * g.dx - eps;
                }
                else // case 2
                {
                    t1 = ((logical_coordinate + 0.5) * g.dx - xp) / x_velocity;
                    s = (logical_coordinate + 0.5) * g.dx + eps;
                }
            }
            else // particle in right half
            {
                if(x_velocity > 0) // case 3
                {
                    t1 = ((logical_coordinate + 1 ) * g.dx - xp)/ x_velocity;
                    s = (logical_coordinate + 1) * g.dx + eps;
                }
                else // case 4
                {
                    t1 = -(xp - (logical_coordinate + 0.5) * g.dx) / x_velocity;
                    s = (logical_coordinate + 0.5) * g.dx - eps;
                }
            }
            // t1 MUST ALWAYS BE POSITIVE ELSE REPOSITIONING LOGIC IS WRONG

            double time_overflow = time_left - t1;
            bool switches_cells = time_overflow  > 0;
            double time_in_this_iteration = switches_cells ? t1 : time_left;
            time_in_this_iteration = (x_velocity == 0.0) ? dt : time_in_this_iteration;

            int logical_coordinate_long = particle_in_left_half ? logical_coordinate: logical_coordinate + 1;
            int logical_coordinate_trans = particle_in_left_half ? logical_coordinate-1: logical_coordinate + 1;

            int sign = (int)(particle_in_left_half) * 2 - 1;
            double distance_to_center = (logical_coordinate + 0.5) * g.dx - xp;
            double s0 = 1 - sign * distance_to_center / g.dx;
            double change_in_coverage = sign * x_velocity * time_in_this_iteration / g.dx;
            double s1 = s0 + change_in_coverage;
            double w = 0.5 * (s0 + s1);

            Array3d j_contribution = v.row(i) * eff_q / dt * time_in_this_iteration;
            double y_contribution_to_current_cell = w * j_contribution(1);
            double z_contribution_to_current_cell = w * j_contribution(2);
            double y_contribution_to_next_cell = (1-w) * j_contribution(1);
            double z_contribution_to_next_cell = (1-w) * j_contribution(2);

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
    g.current_density_x.segment(1, 2) += g.current_density_x.tail(2);
    g.current_density_x.tail(2) = 0;
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

void NonPeriodicSpecies::gather_charge(Grid &g)
{
    gather_charge_computation(g);
}


void NonPeriodicSpecies::gather_current(Grid &g)
{
    gather_current_computation(g);
    g.current_density_yz.topRows(2) = 0;
    g.current_density_yz.bottomRows(2) = 0;
    g.current_density_x(0) = 0;
    g.current_density_x.tail(2) = 0;

}

void Species::apply_particle_bc(Grid &g)
{
    x = x - (g.L * (x/g.L).floor());
}

void NonPeriodicSpecies::apply_particle_bc(Grid &g)
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
        v = new_v;
        E.conservativeResize(N_alive_new, 3);
        B.conservativeResize(N_alive_new, 3);
    }
}


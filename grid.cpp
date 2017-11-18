#include <iostream>
#include <cmath>
#include <eigen3/unsupported/Eigen/FFT>
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

ArrayXd Grid::bincount(ArrayXd cell_numbers, ArrayXd weights, int minlength)
{
    ArrayXd result = ArrayXd::Zero(minlength);
    for (int j = 0; j < cell_numbers.size(); j++)
    {
        int i = cell_numbers[j];
        result[i] += weights[i];
    }
    return result;
}

// CHARGE AND CURRENT DEPOSITION

void Grid::gather_charge(Species s)
{
    ArrayXd logical_coordinates = floor(s.x / dx);
    ArrayXd charge_to_right = (s.x / dx) - logical_coordinates;
    
    ArrayXd charge_hist_to_right = bincount(logical_coordinates+1, charge_to_right, NG+1);
    ArrayXd charge_hist_to_left = bincount(logical_coordinates, 1-charge_to_right, NG+1);
    ArrayXd charge_density(NG);
    charge_density = charge_hist_to_right + charge_hist_to_left;
}

void Grid::gather_charge_periodic(Species s)
{
    gather_charge(s);
    charge_density.head(1) += charge_density.tail(1);
}

typedef Array<bool,Dynamic,1> ArrayXb;

void Grid::gather_current(Species species)
{
    float epsilon = dx * 1e-10;
    for (int i=0; i < species.N; i++)
    {
        float x_velocity = species.v(i,0);
        bool active = species.v.row(i).any();
        float time_left = dt;
        float xp = species.x(i);
        while(active)
        {
            int logical_coordinate = (int)floor(xp/dx);

            bool particle_in_left_half = species.x(i) / dx - logical_coordinate <= 0.5;
            float s, t1;
            if (particle_in_left_half)
            {
                if(x_velocity > 0)
                {
                    t1 = -(xp - logical_coordinate * dx) / x_velocity;
                    s = logical_coordinate * dx - epsilon;
                }
                else
                {
                    t1 = ((logical_coordinate + 0.5) * dx - xp) / x_velocity;
                    s = (logical_coordinate + 0.5) * dx + epsilon;
                }
            }
            else // particle in right half
            {
                if(x_velocity > 0)
                {
                    t1 = ((logical_coordinate + 1 ) * dx)/ x_velocity;
                    s = (logical_coordinate + 1) * dx + epsilon;
                }
                else
                {
                    t1 = -(xp - (logical_coordinate + 0.5) * dx) / x_velocity;
                    s = (logical_coordinate + 0.5) * dx + epsilon;
                }
            }

            float time_overflow = time_left - t1;
            bool switches_cells = time_overflow  > 0;
            float time_in_this_iteration = switches_cells ? t1 : time_left;
            time_in_this_iteration = (x_velocity == 0) ? dt : time_in_this_iteration;

            int logical_coordinate_long = particle_in_left_half ? logical_coordinate: logical_coordinate + 1;
            int logical_coordinate_trans = particle_in_left_half ? logical_coordinate-1: logical_coordinate + 1;

            int sign = (int)(particle_in_left_half) * 2 - 1;
            float distance_to_center = (logical_coordinate + 0.5) * dx - xp;
            float s0 = 1 - sign * distance_to_center / dx;
            float change_in_coverage = sign * x_velocity * time_in_this_iteration / dx;
            float s1 = s0 + change_in_coverage;
            float w = 0.5 * (s0 + s1);

            Array3d j_contribution = species.v.row(i) * species.eff_q / dt * time_in_this_iteration;
            float y_contribution_to_current_cell = w * j_contribution(1);
            float z_contribution_to_current_cell = w * j_contribution(2);
            float y_contribution_to_next_cell = (1-w) * j_contribution(1);
            float z_contribution_to_next_cell = (1-w) * j_contribution(2);

            current_density_x(logical_coordinate_long + 1) += j_contribution(0);
            current_density_yz(logical_coordinate + 2, 0) += y_contribution_to_current_cell;
            current_density_yz(logical_coordinate + 2, 1) += z_contribution_to_current_cell;
            current_density_yz(logical_coordinate_trans + 2, 0) += y_contribution_to_next_cell;
            current_density_yz(logical_coordinate_trans + 2, 1) += z_contribution_to_next_cell;

            active = switches_cells;
            time_left = time_overflow;
            xp = s;
        }
    }
}

void Grid:: gather_current_periodic(Species s)
{
    gather_current(s);
}

// FIELD SOLVERS

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
    // update longitudinal field
    electric_field.col(0) -= dt / epsilon_0 * current_density_x.head(NG+2);

    // update transversal field
    ArrayXd Fplus = 0.5 * (electric_field.col(1) + c * magnetic_field.col(2));
    ArrayXd Fminus = 0.5 * (electric_field.col(1) - c * magnetic_field.col(2));
    ArrayXd Gplus = 0.5 * (electric_field.col(2) + c * magnetic_field.col(1));
    ArrayXd Gminus = 0.5 * (electric_field.col(2) - c * magnetic_field.col(1));

    // propagate forwards
    Fplus.tail(NG-1) = Fplus.head(NG-1) - 0.5 * dt * current_density_yz.block(2, 0, NG-1, 1) / epsilon_0;
    Gplus.tail(NG-1) = Gplus.head(NG-1) - 0.5 * dt * current_density_yz.block(2, 1, NG-1, 1) / epsilon_0;
    // propagate backwards
    Fminus.head(NG-1) = Fminus.tail(NG-1) - 0.5 * dt * current_density_yz.block(2, 0, NG-1, 1) / epsilon_0;
    Gminus.head(NG-1) = Gminus.tail(NG-1) - 0.5 * dt * current_density_yz.block(2, 1, NG-1, 1) / epsilon_0;

    electric_field.col(1) = Fplus + Fminus;
    electric_field.col(2) = Gplus + Gminus;
    magnetic_field.col(1) = (Gplus - Gminus)/c;
    magnetic_field.col(2) = (Fplus - Fminus)/c;
}



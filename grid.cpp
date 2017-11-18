#include <iostream>

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

void Grid::gather_current(Species s)
{
    ArrayXd time_left = ArrayXd::Constant(s.n, dx); // TODO
    float epsilon = dx * 1e-10;
    ArrayXb active = s.v.col(0) != 0; // boolean array
    while(active.any());
    for (int i=0; i < s.N; i++)
    {
        int logical_coordinate = (int)floor(s.x(i), dx);

        // todo switch
        if (s.x(i) / dx - logical_coordinate <= 0.5);// particle_in_left_half:
        {
            float t1 = -(s.x(i)
        }


        //
        
            


    

    return; 
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



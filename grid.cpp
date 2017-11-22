#include <iostream>
#include <cmath>
#include <Eigen/Core>
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

ArrayXd Grid::bincount(const Ref<ArrayXd>& cell_numbers, const Ref<ArrayXd>& weights, int minlength)
{
    ArrayXd result = ArrayXd::Zero(minlength);
    for (int j = 0; j < cell_numbers.size(); j++)
    {
        int i = cell_numbers[j];
        result[i] += weights[i];
    }
    return result;
}

typedef Array<bool,Dynamic,1> ArrayXb;

// FIELD SOLVERS

void Grid::initial_solve(bool neutralize)
{
    size_t dim_x = charge_density.size();
    Eigen::FFT<double> fft;
    Eigen::ArrayXcd rho_F;
    Eigen::ArrayXcd out_final;
    out.setZero(dim_x);
    fft.fwd(rho_F, charge_density);
    ArrayXcd k(dim_x);
    k.setZero(dim_x);
    uint N = floor(dim_x/2) + 1; 
    for (uint i = 0; i < N; i++)
    {
        k.imag()(i) = i;
        k.imag()(i+N) = i - N;
    }

    k.imag() *= 2 * M_PI / (dim_x * dx);
    rho_F /= k;
    fft.inv(rho_F, out_final);
    charge_density = out_final.real();

    /* self.k = 2 * np.pi * fft.fftfreq(self.NG, self.dx) */
    /* self.k[0] = 0.0001 */
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



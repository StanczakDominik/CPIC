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


Grid::Grid(int _NG, float _L, float _c, float _epsilon_0, Temporal &_temporal)
    : c(_c), epsilon_0(_epsilon_0), NG(_NG), L(_L), charge_density(NG+1), current_density_x(NG+3),
    current_density_yz(NG+4, 2), electric_field(NG+2, 3), magnetic_field(NG+2, 3), temporal(_temporal)
{
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
    int dim_x = charge_density.size();
    Eigen::FFT<double> fft;
    Eigen::VectorXcd rho_F(dim_x);
    Eigen::VectorXcd out_final(dim_x);
    rho_F.setZero(dim_x);
    /* cout << rho_F << endl; */
    fft.fwd(rho_F, charge_density.matrix());

    if (neutralize)
    {
        rho_F(0) = 0;
    }

    VectorXcd k(dim_x);
    k.setZero(dim_x);
    int N_fft = floor((dim_x-1)/2) + 1; 
    for (int i = 0; i < N_fft; i++)
    {
        k.imag()(i) = float(i / (float)dim_x);
        if ((i + N_fft < dim_x))
        {
            if ( N_fft % 2 == 0)
            {
                k.imag()(i+N_fft) = float((i - N_fft)/(float)dim_x);
            }
            else
            {
                k.imag()(i+N_fft) = float((i - N_fft + 1)/(float)dim_x);
            }
        }
    }

    k.imag() *= -2 * M_PI / (dim_x * dx);
    /* cout << k << endl; */
    rho_F = rho_F.array() / k.array();
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
    electric_field.col(0) -= temporal.dt / epsilon_0 * current_density_x.head(NG+2);

    // update transversal field
    ArrayXd Fplus = 0.5 * (electric_field.col(1) + c * magnetic_field.col(2));
    ArrayXd Fminus = 0.5 * (electric_field.col(1) - c * magnetic_field.col(2));
    ArrayXd Gplus = 0.5 * (electric_field.col(2) + c * magnetic_field.col(1));
    ArrayXd Gminus = 0.5 * (electric_field.col(2) - c * magnetic_field.col(1));

    // propagate forwards
    Fplus.tail(NG-1) = Fplus.head(NG-1) - 0.5 * temporal.dt * current_density_yz.block(2, 0, NG-1, 1) / epsilon_0;
    Gplus.tail(NG-1) = Gplus.head(NG-1) - 0.5 * temporal.dt * current_density_yz.block(2, 1, NG-1, 1) / epsilon_0;
    // propagate backwards
    Fminus.head(NG-1) = Fminus.tail(NG-1) - 0.5 * temporal.dt * current_density_yz.block(2, 0, NG-1, 1) / epsilon_0;
    Gminus.head(NG-1) = Gminus.tail(NG-1) - 0.5 * temporal.dt * current_density_yz.block(2, 1, NG-1, 1) / epsilon_0;

    electric_field.col(1) = Fplus + Fminus;
    electric_field.col(2) = Gplus + Gminus;
    magnetic_field.col(1) = (Gplus - Gminus)/c;
    magnetic_field.col(2) = (Fplus - Fminus)/c;
}

void Grid::aperiodic_apply_bc(int iteration)
{
    cout << "TODO"; // TODO   
}


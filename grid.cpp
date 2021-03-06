#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <eigen3/unsupported/Eigen/FFT>
#include <Eigen/Dense>
#include "grid.hpp"
#include "H5Cpp.h"
#include <fstream>

using namespace std;
using namespace Eigen;

ArrayXd PythonLinSpaced(int NG, double x_min, double x_max)
{
    // Equivalent to python's numpy.linspace(... endpoint = False;
    ArrayXd x(NG);
    double dx = (x_max-x_min)/NG;
    for (int i = 1; i < NG; i++)
    {
        x(i) = i * dx;
    }
    return x;
}


Grid::Grid(int _NG, double _L, double _c, double _epsilon_0, Temporal &_temporal)
    : c(_c), epsilon_0(_epsilon_0), NG(_NG), L(_L), charge_density(NG+1), current_density_x(NG+3),
    current_density_yz(NG+4, 2), electric_field(NG+2, 3), magnetic_field(NG+2, 3), temporal(_temporal)
{
    x = PythonLinSpaced(NG, 0, L);
    dx = x(1) - x(0);
}

Grid::~Grid()
{
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
        k.imag()(i) = double(i / (double)dim_x);
        if ((i + N_fft < dim_x))
        {
            if ( N_fft % 2 == 0)
            {
                k.imag()(i+N_fft) = double((i - N_fft)/(double)dim_x);
            }
            else
            {
                k.imag()(i+N_fft) = double((i - N_fft + 1)/(double)dim_x);
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

void Grid::apply_bc(double t)
{
    (void)t;
}

NonPeriodicGrid::NonPeriodicGrid(int _NG, double _L, double _c, double _epsilon_0, Temporal &_temporal,
        double _laser_wavelength, double _laser_intensity, double _envelope_center_t, double _envelope_width,
        double _envelope_power)
    : Grid(_NG, _L, _c, _epsilon_0, _temporal),
    laser_omega(2*M_PI*_c/_laser_wavelength),
    laser_amplitude(sqrt(_laser_intensity/(epsilon_0*c))),
    envelope_center_t(_envelope_center_t), envelope_width(_envelope_width),
    envelope_power(_envelope_power),
    _taui(0.5 / pow(log(2), 1./envelope_power) * envelope_center_t),
    _tau(pow(2, 1./envelope_power) * _taui),
    _t_0(_tau * pow(10, 1./envelope_power))
{

}


void NonPeriodicGrid::apply_bc(double t)
{
    // skipping laser phase;
    double wave = laser_amplitude * sin(laser_omega * t); 
    double envelope = exp(-pow((t - _t_0) / _tau, envelope_power));
    double return_value = wave * envelope;
    electric_field(0,1) = return_value;
    magnetic_field(0,2) = return_value / c;
}


string array_filename(string arr, int i)
{
   std::ostringstream stringStream;
   stringStream << arr << "." << i << ".csv";
   return stringStream.str();
}

void Grid::save(int i)
{
    std::ofstream cd(array_filename("charge_density", i));
    cd << charge_density;
    cd.close();
    std::ofstream jx(array_filename("current_density_x", i));
    jx << current_density_x;
    jx.close();
    std::ofstream jyz(array_filename("current_density_yz", i));
    jyz << current_density_yz;
    jyz.close();
    std::ofstream ef(array_filename("electric_field", i));
    ef << electric_field;
    ef.close();
    std::ofstream bf(array_filename("magnetic_field", i));
    bf << magnetic_field;
    bf.close();
}

#ifndef SPECIES_H
#define SPECIES_H
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Grid;

class Species
{
    public:
        float q;
        float m;
        int N;
        int N_alive;
        float scaling;
        float eff_q;
        float eff_m;
        ArrayXd x;
        ArrayX3d v;
        float dt;
        ArrayX3d E;
        ArrayX3d B;
    Species(int N, float q, float m, float scaling);
    double velocity_push();
    void position_push();
    void interpolate_fields(Grid &g);
    void apply_bc(Grid &g);
    void distribute_uniformly(Grid &g, float shift, float start_moat, float end_moat);
    void sinusoidal_position_perturbation(float amplitude, int mode, Grid &g);
    void gather_charge(Grid &g);
    void gather_charge_computation(Grid &g);
    void gather_current(Grid &g);
    void gather_current_computation(Grid &g);
};

class NonPeriodicSpecies : public Species
{
    public:
        using Species::Species;
        void interpolate_fields(Grid &g);
        void apply_bc(Grid &g);
        void gather_charge(Grid &g);
        void gather_current(Grid &g);
};


#endif /* SPECIES_H */

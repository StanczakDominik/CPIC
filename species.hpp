#ifndef SPECIES_H
#define SPECIES_H
#include <Eigen/Dense>
#include "grid.hpp"
using namespace std;
using namespace Eigen;

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
    Species(int N, float q, float m, float scaling, float dt);
    virtual ~Species();
    void position_push();
    double velocity_push(Grid &g);
    virtual void interpolate_fields(Grid &g);
    virtual void apply_particle_bc(Grid &g);
    void distribute_uniformly(Grid &g, float shift, float start_moat, float end_moat);
    void sinusoidal_position_perturbation(float amplitude, int mode, Grid &g);
    virtual void gather_charge(Grid &g);
    void gather_charge_computation(Grid &g);
    virtual void gather_current(Grid &g);
    void gather_current_computation(Grid &g);
};

class NonPeriodicSpecies : public Species
{
    public:
        using Species::Species;
        void apply_particle_bc(Grid &g);
        void interpolate_fields(Grid &g);
        void gather_charge(Grid &g);
        void gather_current(Grid &g);
};


#endif /* SPECIES_H */

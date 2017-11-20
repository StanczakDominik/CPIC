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
    void aperiodic_interpolate_fields(Grid g);
    void periodic_interpolate_fields(Grid g);
    void periodic_apply_bc(Grid g);
    void aperiodic_apply_bc(Grid g);
    void distribute_uniformly(Grid g, float shift, float start_moat, float end_moat);
    void sinusoidal_position_perturbation(float amplitude, int mode, Grid g);
};

#endif /* SPECIES_H */

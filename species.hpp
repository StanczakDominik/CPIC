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
        float scaling;
        float eff_q;
        float eff_m;
        ArrayXd x;
        ArrayX3d v;
        float dt;
    Species(int N, float q, float m, float scaling);
    void velocity_push();
    void position_push();
    ArrayX3d interpolate_electric_field(Grid g);
    ArrayX3d interpolate_magnetic_field(Grid g);
};

#endif /* SPECIES_H */

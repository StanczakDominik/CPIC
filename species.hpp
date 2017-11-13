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
        ArrayX3d E;
        ArrayX3d B;
    Species(int N, float q, float m, float scaling);
    void velocity_push();
    void position_push();
    void interpolate_fields(Grid g);
};

#endif /* SPECIES_H */

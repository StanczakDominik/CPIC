#ifndef SPECIES_H
#define SPECIES_H
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

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
    Species(int N, float q, float m, float scaling);
};

#endif /* SPECIES_H */

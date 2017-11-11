#include <Eigen/Dense>
#include "species.hpp"
#ifndef TEMPORAL_H
#define TEMPORAL_H
using namespace std;
using namespace Eigen;

class Temporal
{
    public:
        float dt;
        float T;
        int NT;

        ArrayXd t;
    Temporal(int _NT, float T);
    /* Temporal(float _dt, float T); */
};

#endif /* TEMPORAL_H */

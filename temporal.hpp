#ifndef TEMPORAL_H
#define TEMPORAL_H
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Temporal
{
    public:
        float dt;
        float T;
        int NT;

        ArrayXd t;
    Temporal(int _NT, float _T);
    Temporal(Temporal &_temporal);
};

#endif /* TEMPORAL_H */

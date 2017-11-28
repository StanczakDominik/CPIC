#ifndef TEMPORAL_H
#define TEMPORAL_H
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Temporal
{
    public:
        double dt;
        double T;
        int NT;

        ArrayXd t;
    Temporal(int _NT, double _T);
    Temporal(double _dt, double _T);
    Temporal(Temporal &_temporal);
};

#endif /* TEMPORAL_H */

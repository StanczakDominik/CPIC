#include <Eigen/Dense>
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
    Temporal(int _NT, float _T);
};

#endif /* TEMPORAL_H */

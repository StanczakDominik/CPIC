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
    Temporal(int _NT, float _T)
    {
        NT = _NT;
        T = _T;
        t = ArrayXd::LinSpaced(0, T, NT);
        dt = t(1) - t(0); // TODO check whether this fits python
    }
    /* Temporal(float _dt, float T); */
};

#endif /* TEMPORAL_H */

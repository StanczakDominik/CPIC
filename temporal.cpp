#include <iostream>
#include <Eigen/Dense>
#include "temporal.hpp"

using namespace std;
using namespace Eigen;
Temporal::Temporal(int _NT, float _T)
    : T(_T), NT(_NT)
{
    t = ArrayXd::LinSpaced(NT, 0, T);
    dt = t(1) - t(0); // TODO check whether this fits python
}

Temporal::Temporal(float _dt, float _T)
    : dt(_dt), T(_T)
{
    NT = (int)(T/dt);
    t = ArrayXd::LinSpaced(NT, 0, T);
}


Temporal::Temporal(Temporal &_temporal)
    : T(_temporal.T), NT(_temporal.NT)
{
    t = _temporal.t;
    dt = _temporal.dt;
}

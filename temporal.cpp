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

Temporal::Temporal(Temporal &_temporal)
{
    Temporal(_temporal.NT, _temporal.T);
}

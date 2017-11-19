#include <iostream>
#include <Eigen/Dense>
#include "temporal.hpp"

using namespace std;
using namespace Eigen;
Temporal::Temporal(int _NT, float _T)
{
    NT = _NT;
    T = _T;
    t = ArrayXd::LinSpaced(0, T, NT);
    dt = t(1) - t(0); // TODO check whether this fits python
}

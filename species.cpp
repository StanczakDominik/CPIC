#include <iostream>
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
    ArrayXd gather_density();
};

Species::Species(int _N, float _q, float _m, float _scaling)
{
    // compute effective charges and masses of macroparticles
    N = _N;
    q = _q;
    m = _m;
    scaling = _scaling;
    eff_q = q * scaling;
    eff_m = m * scaling;
    // allocate position and velocity arrays
    x = ArrayXd(N);
    v = ArrayX3d(N, 3);
}

int main()
{
    Species s = Species(3, 1, 1, 10.0);
    s.x << 1.5, 3.5, 6.5;

    cout << s.x.mean() << endl;
    cout << s.q << endl;
    cout << s.v.size() << endl;
}


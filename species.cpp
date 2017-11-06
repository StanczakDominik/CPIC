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
    Species(int N);
};

Species::Species(int N)
{
    x = ArrayXd(N);
    v = ArrayX3d(N, 3);
}

int main()
{
    Species s = Species(3);
    s.x << 1.5, 3.5, 6.5;
    cout << s.x.mean() << endl;
}


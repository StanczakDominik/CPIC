#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class Grid
{
    public:
        float c;
        float epsilon_0;
        int NG;
        float L;

        ArrayXd x;
        float dx;
        ArrayXd charge_density;
        ArrayXd current_density_x;
        ArrayX2d current_density_yz;
        ArrayX3d electric_field;
        ArrayX3d magnetic_field;
    Grid(int _NG, float _L, float _c, float _epsilon_0);
    ArrayXd gather_density();
};


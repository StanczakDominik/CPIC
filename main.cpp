#include <iostream>
#include "grid.hpp"
#include "species.hpp"
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

void test_grid()
{
    Grid g = Grid(10, 1, 1, 1);

    cout << g.x << endl;
    cout << g.x.mean() << endl;
    cout << g.c << endl;
    cout << g.epsilon_0 << endl;
    cout << g.x.size() << endl;
}

void test_charge_density()
{
    Grid g = Grid(10, 1, 1, 1);
    Species s = Species(10, 1, 1, 1);
    g.gather_charge(s);
    cout << g.charge_density << endl;
}

int main()
{
   Species s = Species(10, 1, 1, 1);
   s.E = ArrayX3d::Zero(s.N, 3);
   s.B = ArrayX3d::Zero(s.N, 3);
   s.B.col(2) = 1.5;
   cout << s.E << endl << "B\n" << s.B << endl;
   cout << "v\n" << s.v << endl << "x\n" << s.x << endl;
   s.velocity_push();
   cout << "v\n" << s.v << endl << "x\n" << s.x << endl;
   
   return 0;
}

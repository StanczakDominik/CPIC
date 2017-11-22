#include <iostream>
#include "grid.hpp"
#include "species.hpp"
#include "simulation.hpp"
#include "temporal.hpp"
/* #include "igl/slice.h" */
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
#include <vector>

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
    s.gather_charge(g);
    cout << g.charge_density << endl;
}

void test_pusher()
{
   Species s = Species(10, 1, 1, 1);
   s.E = ArrayX3d::Zero(s.N, 3);
   s.B = ArrayX3d::Zero(s.N, 3);
   s.B.col(2) = 1.5;
   cout << s.E << endl << "B\n" << s.B << endl;
   cout << "v\n" << s.v << endl << "x\n" << s.x << endl;
   s.velocity_push();
   cout << "v\n" << s.v << endl << "x\n" << s.x << endl;
}

void demonstrate_cross()
{
   Array3d a;  
   Array3d b;
   a <<1,0,0;
   b <<0,1,0;
   cout << a.matrix().cross(b.matrix()) << endl;
}

void demonstrate_cross_rowwise()
{
   ArrayX3d a(4,3);  
   ArrayX3d b(4,3);
   a <<1,0,0,
       0,1,0,
       0,0,1,
       1,0,0;
   b <<0,1,0,
       0,0,1,
       1,0,0,
       0,1,0;
   MatrixX3d am = a.matrix();
   MatrixX3d bm = b.matrix();

   cout << am.rowwise().cross(bm.row(0)) << endl;
}

int main()
{
   cout << "Initializing variables" << endl;
   Species s(1000, 1, 1, 1);
   Grid g(32, 1, 1, 1);
   s.distribute_uniformly(g, 1e-10, 0, 0);
   Temporal temp(1000, 1.0); 
   string file = "filename";
   Simulation sim(temp, g, file, s);
   cout << "Running sim" << endl;
   double runtime = sim.run();
   cout << "Running sim took " << runtime << " seconds" << endl;
   cout << g.charge_density;
}

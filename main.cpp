#include <iostream>
#include "grid.hpp"
#include "species.hpp"
using namespace std;

void test_grid()
{
    Grid g = Grid(10, 1, 1, 1);

    cout << g.x << endl;
    cout << g.x.mean() << endl;
    cout << g.c << endl;
    cout << g.epsilon_0 << endl;
    cout << g.x.size() << endl;
}

int main()
{
    Grid g = Grid(10, 1, 1, 1);
    Species s = Species(10, 1, 1, 1);
    g.gather_charge(s);
    cout << g.charge_density << endl;
    return 0;
}

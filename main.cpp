#include <iostream>
#include "grid.hpp"
#include "species.hpp"
using namespace std;

int main()
{
    Grid g = Grid(10, 1, 1, 1);
    Species s = Species(10, 1, 1, 1);
    g.gather_charge(s);
    cout << g.charge_density << endl;
    return 0;
}

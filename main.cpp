#include <iostream>
#include "grid.hpp"
#include "species.hpp"
using namespace std;

int main()
{
    Grid g = Grid(10, 1, 1, 1);
    Species s = Species(10, 1, 1, 1);
    cout << g.gather_density(s) << endl;
    return 0;
}

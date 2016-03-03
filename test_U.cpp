#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

int main()
{
    ofstream potFile("test_pot.gnu");

    for(double rho=0 ; rho<0.5 ; rho+=0.01)
    {
        potFile << rho << " " << U(rho) << endl;
    }

    return 0;
}

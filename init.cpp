#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"
using namespace std;

int main()
{
    //Initialize random generator
    srand48(time(NULL));

    //Initialize some constants
    double r0 = 1.12; //fm
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    double pf = pow(3./2*M_PI*M_PI*rho0,1./3); //fm-1
    double r = r0 * pow(_A_,1./3);

    //Initialize positions and momenta values
    vector<vector<double> > coords(_NA_, vector<double>(3));
    vector<vector<double> > momenta(_NA_, vector<double>(3));
    sphere(coords, r);
    sphere(momenta, pf);

    //Write the results in in gnu file
    write(coords, "coords.gnu");
    write(momenta, "momenta.gnu");

    return 0;
}

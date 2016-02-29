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
    srand(time(NULL));

    //Initialize some constants
    double r0 = 1.12; //fm
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    double pf = pow(3./2*M_PI*M_PI*rho0,1./3); //fm-1
    unsigned int A = 56;
    unsigned int N = 5;
    unsigned int NA = N*A;
    double r = r0 * pow(A,1./3);

    //Initialize positions and momenta values
    vector<vector<double> > coords(NA, vector<double>(3));
    vector<vector<double> > momenta(NA, vector<double>(3));
    coords = sphere(r, NA);
    momenta = sphere(pf, NA);

    //Write the results in in gnu file
    string const coordsFile("coords.gnu");
    string const momentaFile("momenta.gnu");
    write("coords.gnu", coords);
    write("momenta.gnu", momenta);

    return 0;
}

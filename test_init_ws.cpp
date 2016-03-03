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
    double r = r0 * pow(_A_,1./3);
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    double pf = pow(3./2*M_PI*M_PI*rho0,1./3); //fm-1
    double rms;

    //Initialize positions and momenta values
    vector<vector<double> > coords(_NA_, vector<double>(3));
    vector<vector<double> > momenta(_NA_, vector<double>(3));
    sphere(coords, r);
    sphere(momenta, pf);

    rms = 0;

    for(int i=0 ; i<_NA_ ; i++)
    {
        for(int j=0 ; j<3 ; j++)   
        {
            rms += coords[i][j]*coords[i][j];
        }
    }

    cout << "previous rms = " << sqrt(rms/_NA_) << endl;

    coords_generate(coords, r);
    momenta_generate(coords, momenta);

    rms = 0;

    for(int i=0 ; i<_NA_ ; i++)   
    {
        for(int j=0 ; j<3 ; j++)   
        {
            rms += coords[i][j]*coords[i][j];
        }
    }

    cout << "new rms = " << sqrt(rms/_NA_) << endl;

    return 0;
}

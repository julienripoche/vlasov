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
    double rho0 = 3./4/M_PI/pow(_R0_,3); //fm-3
    double pf = pow(3./2*M_PI*M_PI*rho0,1./3) * _HBAR_C_; //MeV
    double R = _R0_ * pow(_A_,1./3); //fm

    //Initialize positions and momenta values
    vector<vector<double> > coords(_NA_, vector<double>(3));
    vector<vector<double> > momenta(_NA_, vector<double>(3));
    sphere(coords, R);
    sphere(momenta, pf);

    //Write the results in in gnu file
    write(coords, "coords.gnu");
    write(momenta, "momenta.gnu");

    ofstream rpFile("rp.gnu");
    double rMod, pMod;
    for(int i=0 ; i<_NA_ ; i++)
    {
        rMod = coords[i][0]*coords[i][0] + coords[i][1]*coords[i][1] + coords[i][2]*coords[i][2];
        pMod = momenta[i][0]*momenta[i][0] + momenta[i][1]*momenta[i][1] + momenta[i][2]*momenta[i][2];
        rpFile << sqrt(rMod) << " " << sqrt(pMod) << endl;
    }

    return 0;
}

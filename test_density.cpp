#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

//*************************
// RHO MAP NOT FREEZED
//*************************

int main()
{
    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Initialize rho map
    vector<double> rho_map(_BOX_NBR_*_BOX_NBR_*_BOX_NBR_,0);
    rho(rho_map, r);

    //Write density profile in "density.gnu"
    ofstream densityFile("density.gnu");
    for(int i=0 ; i<_BOX_NBR_ ; i++)
    {
        for(int j=0 ; j<_BOX_NBR_ ; j++)
        {
            densityFile << i*_L0_ << " " << j*_L0_ << " " << rho_map[key(i,j,_BOX_NBR_/2,_BOX_NBR_)] << endl;
        }
        densityFile << endl;
    }

    //Write potential profile in "ubar.gnu" 
    ofstream ubarFile("ubar.gnu");
    vector<double> ubar_map(_BOX_NBR_*_BOX_NBR_*_BOX_NBR_,0);
    vector<double> r0(3);
    for(int i=0 ; i<_BOX_NBR_ ; i++)
    {
        r0[0] = (i-_BOX_NBR_/2.)*_L0_;
        for(int j=0 ; j<_BOX_NBR_ ; j++)
        {
            r0[1] = (j-_BOX_NBR_/2.)*_L0_;
            
            for(int k=0 ; k<_BOX_NBR_ ; k++)
            {
                r0[2] = (k-_BOX_NBR_/2.)*_L0_;
                ubar_map[key(i,j,k,_BOX_NBR_)] = get_ubar(rho_map, r0);
            }
            ubarFile << i*_L0_ << " " << j*_L0_ << " " << ubar_map[key(i,j,_BOX_NBR_/2,_BOX_NBR_)] << endl;
        }
        ubarFile << endl;
    }

    return 0;
}

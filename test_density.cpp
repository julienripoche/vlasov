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
    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Useful variables
    int n_box = 30;
    double L = 20.;
    double l0 = L/n_box;

    //Initialize rho map
    vector<double> rho_map(n_box*n_box*n_box,0);
    rho(rho_map, r, l0, n_box);

    //Write density profile in "density.gnu" 
    ofstream densityFile("density.gnu");
    for(int i=0 ; i<n_box ; i++)
    {
        for(int j=0 ; j<n_box ; j++)
        {
            densityFile << i << " " << j << " " << rho_map[key(i,j,n_box/2,n_box)] << endl;
        }
    }

    return 0;
}

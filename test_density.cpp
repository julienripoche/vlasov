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
    read("coords.gnu", r);
    read("momenta.gnu", p);

    //Initialize rho map
    int n = 30;
    vector<double> rho(n*n*n,0);
    rho_map(rho, r, 1., 20, 30);

    //Write density profile in "density.gnu" 
    ofstream densityFile("density.gnu");
    for(int i=0 ; i<n ; i++)
    {
        for(int j=0 ; j<n ; j++)
        {
            densityFile << i << " " << j << " " << rho[key(i,j,n/2,n)] << endl;
        }
    }

    return 0;
}

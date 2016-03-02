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

    //Usefel variables
    int n_box = 30;
    int sigma_nbr = 2;
    double L = 15;
    double l0 = L/n_box;

    //Initialize rho map
    vector<double> rho_map(n_box*n_box*n_box,0);
    rho(rho_map, r, l0, n_box);

    //Initialize ubar map
    //Write density profile in "ubar.gnu" 
    ofstream ubarFile("ubar.gnu");
    vector<double> ubar_map(n_box*n_box*n_box,0);
    vector<double> r0(3);
    for(int i=0 ; i<n_box ; i++)
    {
        r0[0] = (i-n_box/2.)*l0;
        for(int j=0 ; j<n_box ; j++)
        {
            r0[1] = (j-n_box/2.)*l0;
            for(int k=0 ; k<n_box ; k++)
            {
                r0[2] = (k-n_box/2.)*l0;
                ubar_map[key(i,j,k,n_box)] = get_ubar(rho_map, r0, l0, n_box, sigma_nbr);
            }
            ubarFile << i*l0 << " " << j*l0 << " " << ubar_map[key(i,j,n_box/2,n_box)] << endl;
        }
        ubarFile << endl;
    }

    return 0;
}

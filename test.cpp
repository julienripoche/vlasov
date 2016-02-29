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

    unsigned int n = 30;

    //Initialize strength
    unsigned int NA = r.size();
    vector<vector<double> > F(NA, vector<double>(3,0)); 
    vector<double> rho(n*n*n,0);

    rho_map(rho, r, 1., 20, 30, 5);

    //for(unsigned int i=0 ; i<NA ; i++)
    //{
        F[0] = minus_gradU(r[0], rho, 2, n, 20, 0.4);
    //}

    cout << F[0][0] << " " << F[0][1] << " " << F[0][2] << endl;

    return 0;
}

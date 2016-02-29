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
    //Initialize random generator
    srand(time(NULL));

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

    rho_map(rho, r, 1., 20, 30);

    double sum;
 
    ofstream densityFile("density.gnu");
    for(unsigned int i=0 ; i<n ; i++)
    {
        for(unsigned int j=0 ; j<n ; j++)
        {
            densityFile << i << " " << j << " " << rho[key(i,j,n/2,n)] << endl;
        }
    }

    cout << "aaa" << endl;

    //for(unsigned int i=0 ; i<NA ; i++)
    //{
    cout << r[1119][2] << endl;

    vector<double> truc(3,0);
    F[0] = truc;

    cout << "0 value" << endl;

    F[0] = minus_gradU2(r[0], rho, 2, n, 20, 0.4);
    //}

    cout << "bbb" << endl;

    return 0;
}

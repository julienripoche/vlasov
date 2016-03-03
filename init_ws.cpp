#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"
#include "wsRadius.h"

using namespace std;

int main()
{
    //Initialize random generator
    srand(time(NULL));

    //Initialize ws radius
    //double R = _R0_ * pow(_A_,1./3);
    WsRadiusFinder WsRadius(_SIGMA_/2);
    double R = WsRadius.run();

    //Initialize positions and momenta values
    vector<vector<double> > coords(_NA_, vector<double>(3));
    vector<vector<double> > momenta(_NA_, vector<double>(3));
    coords_generate(coords, R);
    momenta_generate(coords, momenta);

    //Write the results in in gnu file
    write(coords, "coords.gnu");
    write(momenta, "momenta.gnu");

    //Check the relation between r and p
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

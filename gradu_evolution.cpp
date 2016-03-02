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

    //Initialize strength
    vector<vector<double> > F(_NA_, vector<double>(3,0));
    for(int i=0 ; i<_NA_ ; i++)
    {
        minus_gradU(F[i], rho_map, r[i]);
    }

    // Initialize useful variables
    int n_ite = 100;
    double r_modulus;
    double p_modulus;
    double r_rms;
    double p_rms;

    //Open file to write results
    ofstream partFile("particle.gnu");
    ofstream rmsFile("rms.gnu");

    for(int i=0 ; i<n_ite ; i++)
    {
        //Initialize r and p rms
        r_rms = 0;
        p_rms = 0;

        //Loop over all test particles
        for(int j=0 ; j<_NA_ ; j++)
        {
            //Initialize r and p modulus value
            r_modulus = 0;
            p_modulus = 0;

            //Loop over cartesian coordinates
            for(int k=0 ; k<3 ; k++)
            {
                //Differential system resolution
                p[j][k] += 1./2*F[j][k]*_DT_;
                r[j][k] += 1./_M_*p[j][k]*_DT_;
            }
        }

        rho(rho_map, r);

        for(int j=0 ; j<_NA_ ; j++)
        {
            minus_gradU(F[j], rho_map, r[j]);

            //Loop over cartesian coordinates
            for(int k=0 ; k<3 ; k++)
            { 
                p[j][k] += 1./2*F[j][k]*_DT_;

                //Rms
                r_rms += r[j][k]*r[j][k];
                p_rms += p[j][k]*p[j][k];

                //For one test particle
                if(j==1)
                {
                    r_modulus += r[j][k]*r[j][k];
                    p_modulus += p[j][k]*p[j][k];
                }
            }

            //Write modulus of one test particle
            if(j==1)
            {
                partFile << i*_DT_ << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        //To follow the progession of the evolution
        cout << i << "/" << n_ite << endl;

        //Write rms values in gnu files
        rmsFile << i*_DT_ << " " << sqrt(r_rms/_NA_) << " " << sqrt(p_rms/_NA_) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords.gnu");
    write(p, "after_momenta.gnu");

    return 0;
}

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
    double omega = 0.1; //fm-1

    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Initialize strength
    vector<vector<double> > F(_NA_, vector<double>(3,0));
    for(int i=0 ; i<_NA_ ; i++)
    {
        for(int j=0 ; j<3 ; j++)
        {
            F[i][j] = -_M_*omega*omega*r[i][j];
        }
    }

    // Initialize useful variables
    int n_ite = 1000;
    double r_modulus;
    double p_modulus;
    double r_rms;
    double p_rms;

    //Open file to write results
    ofstream partFile("particle_harm.gnu");
    ofstream rmsFile("rms_harm.gnu");

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
                F[j][k] = -_M_*omega*omega*r[j][k]; //Using harmonic oscillator
                p[j][k] += 1./2*F[j][k]*_DT_;

                //Rms and total momentum calculation
                r_rms += r[j][k]*r[j][k];
                p_rms += p[j][k]*p[j][k];

                //For one test particle
                if(j==2)
                {
                    r_modulus += r[j][k]*r[j][k];
                    p_modulus += p[j][k]*p[j][k];
                }
            }

            //Write modulus of one test particle
            if(j==2)
            {
                partFile << i*_DT_ << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        //Write r and p rms
        rmsFile << i*_DT_ << " " << sqrt(r_rms/_NA_) << " " << sqrt(p_rms/_NA_) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords_harm.gnu");
    write(p, "after_momenta_harm.gnu");

    return 0;
}

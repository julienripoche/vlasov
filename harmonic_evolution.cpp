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
    //Initialize some constants
    double m_proton = 938.3; //MeV
    double m_neutron = 939.6; //MeV
    double hbar_c = 197.3; //MeV.fm
    double m = (m_proton+m_neutron)/2/hbar_c; //fm-1
    double Dt = 0.01; //fm
    double omega = 1.; //fm-1

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
            F[i][j] = -m*omega*omega*r[i][j];
        }
    }

    // Initialize useful variables
    int n_ite = 1000;
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
                p[j][k] += 1./2*F[j][k]*Dt;
                r[j][k] += 1./m*p[j][k]*Dt;
                F[j][k] = -m*omega*omega*r[j][k]; //Using harmonic oscillator
                p[j][k] += 1./2*F[j][k]*Dt;

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
                partFile << i*Dt << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        //Write r and p rms
        rmsFile << i*Dt << " " << sqrt(r_rms/_NA_) << " " << sqrt(p_rms/_NA_) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords.gnu");
    write(p, "after_momenta.gnu");

    return 0;
}

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
    read("coords.gnu", r);
    read("momenta.gnu", p);

    //Initialize strength
    unsigned int NA = r.size();
    vector<vector<double> > F(NA, vector<double>(3,0));
    for(unsigned int i=0 ; i<NA ; i++)
    {
        for(unsigned int j=0 ; j<3 ; j++)
        {
            F[i][j] = -m*omega*omega*r[i][j];
        }
    }

    // Initialize useful variables
    unsigned int n_ite = 1000;
    double r_modulus;
    double p_modulus;
    double rms;
    double total_momentum[3];

    //Open file to write results
    ofstream particleFile("one_particle_oscillation.gnu");
    ofstream rmsFile("rms.gnu");
    ofstream momentumConsFile("momentum_cons.gnu");

    for(unsigned int i=0 ; i<n_ite ; i++)
    {
        //Initialize RMS and total momentum
        rms = 0;
        for(unsigned int j=0 ; j<3 ; j++)
        {
            total_momentum[j] = 0;
        }

        //Loop over all test particles
        for(unsigned int j=0 ; j<NA ; j++)
        {
            //Initialize r and p modulus value
            r_modulus = 0;
            p_modulus = 0;

            //Loop over cartesian coordinates
            for(unsigned int k=0 ; k<3 ; k++)
            {
                //Differential system resolution
                p[j][k] += 1./2*F[j][k]*Dt;
                r[j][k] += 1./m*p[j][k]*Dt;
                F[j][k] = -m*omega*omega*r[j][k]; //Using harmonic oscillator
                p[j][k] += 1./2*F[j][k]*Dt;

                //Rms and total momentum calculation
                rms += r[j][k]*r[j][k];
                total_momentum[k] += p[j][k];

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
                particleFile << i*Dt << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        //Write rms and momenta values in gnu files
        rmsFile << sqrt(rms/NA) << endl;
        momentumConsFile << total_momentum[0] << " " << total_momentum[1] << " " << total_momentum[2] << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write("after_coords.gnu", r);
    write("after_momenta.gnu", p);

    return 0;
}

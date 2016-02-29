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
    double omega = 0.5; //fm-1

    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read("coords.gnu", r);
    read("momenta.gnu", p);

    //Initialize strength
    unsigned int part_nbr = r.size();
    vector<vector<double> > F(part_nbr, vector<double>(3,0));
    for(unsigned int i=0 ; i<part_nbr ; i++)
    {
        for(unsigned int j=0 ; j<3 ; j++)
        {
            F[i][j] = -m*omega*omega*r[i][j];
        }
    }

    // Initialize useful variables
    unsigned int n_ite = 10000;
    double r_modulus;
    double rms;
    double total_momentum[3];
    double Dp_half;

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
        for(unsigned int j=0 ; j<part_nbr ; j++)
        {
            //Initialize r modulus value
            r_modulus = 0;

            //Loop over cartesian coordinates
            for(unsigned int k=0 ; k<3 ; k++)
            {
                //Differential system resolution
                Dp_half = 1./2*F[j][k]*Dt;
                p[j][k] += Dp_half;
                r[j][k] += 1./m*p[j][k]*Dt;
                F[j][k] = -m*omega*omega*r[j][k]; //Using harmonic oscillator
                p[j][k] += Dp_half;

                //Rms and total momentum calculation
                rms += r[j][k]*r[j][k];
                total_momentum[k] += p[j][k];

                //For the first test particle
                if(j==0)
                {
                    r_modulus += r[j][k]*r[j][k];
                }
            }

            //Write modulus of the first test particle
            if(j==0)
            {
                particleFile << i*Dt << " " << sqrt(r_modulus) << endl;
            }
        }

        //Write rms and momenta values in gnu files
        rmsFile << sqrt(rms)/part_nbr << endl;
        momentumConsFile << total_momentum[0] << " " << total_momentum[1] << " " << total_momentum[2] << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write("after_coords.gnu", r);
    write("after_momenta.gnu", p);

    return 0;
}

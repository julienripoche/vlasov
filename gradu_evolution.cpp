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
    double Dt = 0.005; //fm

    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read("coords.gnu", r);
    read("momenta.gnu", p);

    //Initialize strength
    int NA = r.size();
    vector<vector<double> > F(NA, vector<double>(3,0)); 

    //Initialize rho map
    int n = 30;
    vector<double> rho(n*n*n,0);
    rho_map(rho, r, 1., 20, n, 5);

    for(int i=0 ; i<NA ; i++)
    {
        F[i] = minus_gradU(r[i], rho, 2, n, 20, 1.);
    }

    // Initialize useful variables
    int n_ite = 100;
    double r_modulus;
    double p_modulus;
    double r_rms;
    double p_rms;

    //Open file to write results
    ofstream particleFile("one_particle_oscillation.gnu");
    ofstream rmsFile("rms.gnu");

    for(int i=0 ; i<n_ite ; i++)
    {
        //Initialize rms
        r_rms = 0;
        p_rms = 0;

        //Loop over all test particles
        for(int j=0 ; j<NA ; j++)
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
            }
        }

        rho_map(rho, r, 1., 20, n, 5);

        for(int j=0 ; j<NA ; j++)
        {
            F[j] = minus_gradU(r[j], rho, 2, n, 20, 1.);

            //Loop over cartesian coordinates
            for(int k=0 ; k<3 ; k++)
            { 
                p[j][k] += 1./2*F[j][k]*Dt;

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
                particleFile << i*Dt << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        cout << endl <<i << "/" << n_ite << endl;

        //Write rms values in gnu files
        rmsFile << i*Dt << " " << sqrt(r_rms/NA) << " " << sqrt(p_rms/NA) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write("after_coords.gnu", r);
    write("after_momenta.gnu", p);

    return 0;
}

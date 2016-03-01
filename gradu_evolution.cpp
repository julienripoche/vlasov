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
    double Dt = 0.05; //fm

    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Usefel variables
    int n_box = 20;
    int nbr_sigma = 2;
    double L = 15;
    double l0 = L/n_box;

    //Initialize rho map
    vector<double> rho_map(n_box*n_box*n_box,0);
    rho(rho_map, r, l0, n_box);

    //Initialize strength
    vector<vector<double> > F(_NA_, vector<double>(3,0));
    for(int i=0 ; i<_NA_ ; i++)
    {
        minus_gradU(F[i], rho_map, r[i], nbr_sigma, n_box, l0);
    }

    // Initialize useful variables
    int n_ite = 200;
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
            }
        }

        rho(rho_map, r, l0, n_box);

        for(int j=0 ; j<_NA_ ; j++)
        {
            minus_gradU(F[j], rho_map, r[j], nbr_sigma, n_box, l0);

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
                partFile << i*Dt << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << endl;
            }
        }

        //To follow the progession of the evolution
        cout << i << "/" << n_ite << endl;

        //Write rms values in gnu files
        rmsFile << i*Dt << " " << sqrt(r_rms/_NA_) << " " << sqrt(p_rms/_NA_) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords.gnu");
    write(p, "after_momenta.gnu");

    return 0;
}

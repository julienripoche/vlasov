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
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Usefel variables
    int N = 5;
    int n_box = 30;
    int nbr_sigma = 2;
    double sigma = 1;
    double L = 20;
    double l0 = L/n_box;

    //Initialize rho map
    vector<double> rho(n_box*n_box*n_box,0);
    rho_map(rho, r, sigma, l0, n_box, N);

    //Initialize strength
    int NA = r.size();
    vector<vector<double> > F(NA, vector<double>(3,0));
    for(int i=0 ; i<NA ; i++)
    {
        minus_gradU(F[i], r[i], rho, nbr_sigma, n_box, l0, sigma);
    }

    // Initialize useful variables
    int n_ite = 5;
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

        rho_map(rho, r, sigma, l0, n_box, N);

        for(int j=0 ; j<NA ; j++)
        {
            minus_gradU(F[j], r[j], rho, nbr_sigma, n_box, l0, sigma);

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
        rmsFile << i*Dt << " " << sqrt(r_rms/NA) << " " << sqrt(p_rms/NA) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords.gnu");
    write(p, "after_momenta.gnu");

    return 0;
}

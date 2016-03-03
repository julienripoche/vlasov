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
    //Initialize positions and momenta values for particle 1
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Initialize positions and momenta values for particle 2
    vector<vector<double> > r2;
    vector<vector<double> > p2;
    read(r2, "coords2.gnu");
    read(p2, "momenta2.gnu");

    //Add particle 2 values to particle 1 vectors
    int size = r2.size();
    for(int i=0 ; i<size ; i++)
    {
        r.push_back(r2[i]);
        p.push_back(p2[i]);
    }
    
    cout << "copy " << r.size() << " " << p.size() << endl;

    return 0;

    cout << "elsewhere" << endl;

    //Initialize rho map
    vector<double> rho_map(_BOX_NBR_X_*_BOX_NBR_Y_*_BOX_NBR_Z_,0);
    //rho2(rho_map, r);

    double px = 0;
    double py = 0;
    double pz = 0;
    double ptot = 0;

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

    //Write potential profile in "ubar.gnu" 
    ofstream ubarFile("ubar.gnu");
    vector<double> ubar_map(_BOX_NBR_*_BOX_NBR_*_BOX_NBR_,0);
    vector<double> r0(3);
    for(int i=0 ; i<_BOX_NBR_ ; i++)
    {
        r0[0] = (i-_BOX_NBR_/2.)*_L0_;
        for(int j=0 ; j<_BOX_NBR_ ; j++)
        {
            r0[1] = (j-_BOX_NBR_/2.)*_L0_;
            
            for(int k=0 ; k<_BOX_NBR_ ; k++)
            {
                r0[2] = (k-_BOX_NBR_/2.)*_L0_;
                ubar_map[key(i,j,k,_BOX_NBR_)] = get_ubar(rho_map, r0);
            }
            ubarFile << i*_L0_ << " " << j*_L0_ << " " << ubar_map[key(i,j,_BOX_NBR_/2,_BOX_NBR_)] << endl;
        }
        ubarFile << endl;
    }

    //Initialize strength
    vector<vector<double> > F(_NA_, vector<double>(3,0));
    for(int i=0 ; i<_NA_ ; i++)
    {
        minus_gradU(F[i], rho_map, r[i]);
    }

    // Initialize useful variables
    int n_ite = 1000;
    double r_modulus;
    double p_modulus;
    double gradu_modulus;
    double r_rms;
    double p_rms;
    double gradu_rms;
    char densityFileName[50];
    char ubarFileName[50];

    //Open file to write results
    ofstream partFile("particle.gnu");
    ofstream rmsFile("rms.gnu");

    for(int i=0 ; i<n_ite ; i++)
    {
        //Initialize r and p rms
        r_rms = 0;
        p_rms = 0;
        gradu_rms = 0;

        //Loop over all test particles
        for(int j=0 ; j<_NA_ ; j++)
        {
            //Initialize r and p modulus value
            r_modulus = 0;
            p_modulus = 0;
            gradu_modulus = 0;

            //Loop over cartesian coordinates
            for(int k=0 ; k<3 ; k++)
            {
                //Differential system resolution
                p[j][k] += 1./2*F[j][k]*_DT_;
                r[j][k] += 1./_M_*p[j][k]*_DT_;
            }
        }

        rho(rho_map, r);

        //Write density profile
        sprintf(densityFileName, "density/density%d.gnu", i);
        ofstream densityFile(densityFileName);
        for(int i2=0 ; i2<_BOX_NBR_ ; i2++)
        {
            for(int j2=0 ; j2<_BOX_NBR_ ; j2++)
            {
                densityFile << i2*_L0_ << " " << j2*_L0_ << " " << rho_map[key(i2,j2,_BOX_NBR_/2,_BOX_NBR_)] << endl;
            }
            densityFile << endl;
        }

        /*
        //Write potential profile
        sprintf(ubarFileName, "ubar/ubar%d.gnu", i);
        ofstream ubarFile(ubarFileName);
        vector<double> ubar_map(_BOX_NBR_*_BOX_NBR_*_BOX_NBR_,0);
        vector<double> r0(3);
        for(int i2=0 ; i2<_BOX_NBR_ ; i2++)
        {
            r0[0] = (i2-_BOX_NBR_/2.)*_L0_;
            for(int j2=0 ; j2<_BOX_NBR_ ; j2++)
            {
                r0[1] = (j2-_BOX_NBR_/2.)*_L0_;

                for(int k2=0 ; k2<_BOX_NBR_ ; k2++)
                {
                    r0[2] = (k2-_BOX_NBR_/2.)*_L0_;
                    ubar_map[key(i2,j2,k2,_BOX_NBR_)] = get_ubar(rho_map, r0);
                }
                ubarFile << i2*_L0_ << " " << j2*_L0_ << " " << ubar_map[key(i2,j2,_BOX_NBR_/2,_BOX_NBR_)] << endl;
            }
            ubarFile << endl;
        }
        */

        px = 0;
        py = 0;
        pz = 0;
        ptot = 0;

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
                gradu_rms += F[j][k]*F[j][k];

                //For one test particle
                if(j==1)
                {
                    r_modulus += r[j][k]*r[j][k];
                    p_modulus += p[j][k]*p[j][k];
                    gradu_modulus += F[j][k]*F[j][k];
                }
            }

            px += p[i][0];
            py += p[i][1];
            pz += p[i][2];

            //For one test particle
            if(j==1 && i%5==0)
            {
                //cout << "x = " << r[j][0] << "  y = " << r[j][1] << "  z = " << r[j][2] << endl;
                //cout << "px = " << p[j][0] << " py = " << p[j][1] << " pz = " << p[j][2] << endl;
                //cout << "Fx = " << F[j][0] << " Fy = " << F[j][1] << " Fz = " << F[j][2] << endl;
            }

            //Write modulus of one test particle
            if(j==1)
            {
                partFile << i*_DT_ << " " << sqrt(r_modulus) << " " << sqrt(p_modulus) << " " << sqrt(gradu_modulus) << endl;
            }
        }

        ptot = sqrt(px*px + py*py + pz*pz)/_NA_;
        //cout << "px = " << px/_NA_ << " py = " << py/_NA_ << " pz = " << pz/_NA_ << " ptot = " << ptot <<  endl;

        //To follow the progession of the evolution
        cout << i << "/" << n_ite << endl;

        //Write rms values in gnu files
        rmsFile << i*_DT_ << " " << sqrt(r_rms/_NA_) << " " << sqrt(p_rms/_NA_) << " " << sqrt(gradu_rms/_NA_) << endl;
    }

    //Write the coordinates and momenta results in gnu files
    write(r, "after_coords.gnu");
    write(p, "after_momenta.gnu");

    return 0;
}

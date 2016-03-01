#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

double alea()
{
    return (double)rand()/RAND_MAX;
}

void sphere(vector<vector<double> > &r, double radius_max)
{
    double radius, theta, phi;
    int nbr_points = r.size();

    for(int i=0 ; i<nbr_points ; i++)
    {
        radius = radius_max*pow(alea(),1./3);
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        r[i][0] = radius*sin(theta)*cos(phi);
        r[i][1] = radius*sin(theta)*sin(phi);
        r[i][2] = radius*cos(theta);
    }
}

void write(vector<vector<double> > &data, string const nomFichier)
{
    ofstream monFlux(nomFichier.c_str());
    int data_size = data.size();

    if(monFlux)
    {
        for(int i=0 ; i<data_size ; i++)
        {
            monFlux << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;
        }
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en écriture." << endl;
    }
}

void read(vector<vector<double> > &data, string const nomFichier)
{
    ifstream monFlux(nomFichier.c_str());
    vector<double> r(3);

    while(monFlux >> r[0] >> r[1] >> r[2])
    {
        data.push_back(r);
    }
}

int key(int x, int y, int z, int N)
{
    return x*N*N + y*N + z;
}

double module(vector<double> &r_real, vector<double> &r_box)
{
    double sqr_sum = 0;
    for(int i=0 ; i<3 ; i++)
    {
        sqr_sum += (r_real[i] - r_box[i]) * (r_real[i] - r_box[i]);
    }
    return sqrt(sqr_sum);
}

double gaussian(vector<double> &r_real, vector<double> &r_box)
{
    double r = module(r_real, r_box);
    return 1/pow(sqrt(2*M_PI)*_SIGMA_,3)*exp(-r*r/2/_SIGMA_/_SIGMA_);
}

double U(double rho)
{
    double r0 = 1.12; //fm
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    return -356*rho/rho0 + 303*pow(rho/rho0,7./6);
}

void rho(vector<double> &rho_map, vector<vector<double> > &coords, double l0, int box_size)
{
    //Initialize some variables
    vector<double> r(3);

    //Loop over the grid
    for(int x=0 ; x<box_size ; x++)
    {
        r[0] = (x-box_size/2.)*l0;
        for(int y=0 ; y<box_size ; y++)
        {
            r[1] = (y-box_size/2.)*l0;
            for(int z=0 ; z<box_size ; z++)
            {
                r[2] = (z-box_size/2.)*l0;
                for(int i=0 ; i<_NA_ ; i++)
                {
                    rho_map[key(x,y,z,box_size)] += gaussian(r, coords[i]);
                }
                rho_map[key(x,y,z,box_size)] /= _N_;
            }
        }
    }
}

void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r, int nbr_sigma, int N, double l0)
{
    //Define some useful variables
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    int nbr_cells = floor(nbr_sigma*_SIGMA_/l0);

    int x1, y1, z1;

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        //box_init[i] = (floor(r[i]/l0) + 0.5)*l0;
        box_init[i] = floor(r[i]/l0);
        gradu[i] = 0;
    }

    //Loop over considered cells
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        box[0] = (box_init[0] + x + 0.5)*l0;
        for(int y=-nbr_cells ; y<=nbr_cells ; y++)
        {
            box[1] = (box_init[1] + y + 0.5)*l0;
            for(int z=-nbr_cells ; z<=nbr_cells ; z++)
            {
                box[2] = (box_init[2] + z + 0.5)*l0;
                
                gaus = gaussian(r, box);

                x1 = (box_init[0] + N/2 + x)%N;
                y1 = (box_init[1] + N/2 + y)%N;
                z1 = (box_init[2] + N/2 + z)%N;

                rho = rho_map[key(x1,y1,z1,N)];

                gradu[0] += (r[0] - (x1 + 0.5)*l0) * gaus * U(rho);
                gradu[1] += (r[1] - (y1 + 0.5)*l0) * gaus * U(rho);
                gradu[2] += (r[2] - (z1 + 0.5)*l0) * gaus * U(rho);

                /*
                for(int i=0 ; i<3 ; i++)
                {
                    gradu[i] += (r[i] - (box_init[i]+[i]+0.5)*l0) * gaus * U(rho);
                }*/
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[i] *= l0*l0*l0/_SIGMA_/_SIGMA_;
    }
}


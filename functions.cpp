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
    return double(rand())/RAND_MAX;
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

double rho_ws(double r)
{
    double r0 = 1.12;
    double R = r0 * pow(_A_, 1./3);
    double rho0 = 3./4/M_PI/pow(r0,3);
    return rho0/(1+exp((r-R)/_SIGMA_));
}

void coords_generate(vector<vector<double> > &r, double radius_max)
{
    double r0 = 1.12;
    double rho0 = 3./4/M_PI/pow(r0,3);
    double radius, theta, phi;
    int nbr_points = r.size();

    for(int i=0 ; i<nbr_points ; i++)
    {
        do
        {
            radius = radius_max*pow(alea(),1./3);
        }
        while(rho_ws(radius)/rho0 < alea());
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        r[i][0] = radius*sin(theta)*cos(phi);
        r[i][1] = radius*sin(theta)*sin(phi);
        r[i][2] = radius*cos(theta);
    }
}

double fermi_momentum(double rho)
{
    return pow(3./2*M_PI*M_PI*rho,1./3);
}

void momenta_generate(vector<vector<double> > &r, vector<vector<double> > &p)
{
    double radius, theta, phi;
    int nbr_points = r.size();
    double pf;
    double module;

    for(int i=0 ; i<nbr_points ; i++)
    {
        module = sqrt(r[i][0]*r[i][0] + r[i][1]*r[i][1] + r[i][2]*r[i][2]);
        pf = fermi_momentum(rho_ws(module));
        radius = pf*pow(alea(),1./3);
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        p[i][0] = radius*sin(theta)*cos(phi);
        p[i][1] = radius*sin(theta)*sin(phi);
        p[i][2] = radius*cos(theta);
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

void rho(vector<double> &rho_map, vector<vector<double> > &coords, double l0, int box_nbr)
{
    //Initialize some variables
    vector<double> r(3);

    //Loop over the grid
    for(int x=0 ; x<box_nbr ; x++)
    {
        r[0] = (x-box_nbr/2.)*l0;
        for(int y=0 ; y<box_nbr ; y++)
        {
            r[1] = (y-box_nbr/2.)*l0;
            for(int z=0 ; z<box_nbr ; z++)
            {
                r[2] = (z-box_nbr/2.)*l0;
                for(int i=0 ; i<_NA_ ; i++)
                {
                    rho_map[key(x,y,z,box_nbr)] += gaussian(r, coords[i]);
                }
                rho_map[key(x,y,z,box_nbr)] /= _N_;
            }
        }
    }
}

void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r, double l0, int box_nbr, int sigma_nbr)
{
    //Define some useful variables
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    int nbr_cells = floor(sigma_nbr*_SIGMA_/l0);

    int x1, y1, z1;

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
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

                x1 = (box_init[0] + x + box_nbr/2 + box_nbr)%box_nbr;
                y1 = (box_init[1] + y + box_nbr/2 + box_nbr)%box_nbr;
                z1 = (box_init[2] + z + box_nbr/2 + box_nbr)%box_nbr;

                rho = rho_map[key(x1,y1,z1,box_nbr)];

                gradu[0] += (r[0] - (box_init[0] + x + 0.5)*l0) * gaus * U(rho);
                gradu[1] += (r[1] - (box_init[1] + y + 0.5)*l0) * gaus * U(rho);
                gradu[2] += (r[2] - (box_init[2] + z + 0.5)*l0) * gaus * U(rho);
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[i] *= l0*l0*l0/_SIGMA_/_SIGMA_;
    }
}

double get_ubar(vector<double> &rho_map, vector<double> &r, double l0, int box_nbr, int sigma_nbr)
{
    //Define some useful variables
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    int nbr_cells = floor(sigma_nbr*_SIGMA_/l0);
    double ubar;
    int x1, y1, z1;

    //Initialize box coordinates
    ubar = 0;
    for(int i=0 ; i<3 ; i++)
    {
        box_init[i] = floor(r[i]/l0);
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

                //(... + box_nbr)%box_nbr; to avoid negative number
                x1 = (box_init[0] + x + box_nbr/2 + box_nbr)%box_nbr;
                y1 = (box_init[1] + y + box_nbr/2 + box_nbr)%box_nbr;
                z1 = (box_init[2] + z + box_nbr/2 + box_nbr)%box_nbr;
                rho = rho_map[key(x1,y1,z1,box_nbr)];

                ubar += gaus * U(rho);
            }
        }
    }

    return ubar*l0*l0*l0;
}



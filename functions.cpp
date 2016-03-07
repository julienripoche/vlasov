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
    //double radius, theta, phi;
    int nbr_points = r.size();

    /*for(int i=0 ; i<nbr_points ; i++)
    {
        radius = radius_max*pow(alea(),1./3);
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        r[i][0] = radius*sin(theta)*cos(phi);
        r[i][1] = radius*sin(theta)*sin(phi);
        r[i][2] = radius*cos(theta);
    }*/

    for(int i=0 ; i<nbr_points ; i++)
    {
        do
        {
        r[i][0] = radius_max*(1-2*alea());
        r[i][1] = radius_max*(1-2*alea());
        r[i][2] = radius_max*(1-2*alea());
        } while(sqrt(r[i][0]*r[i][0]+r[i][1]*r[i][1]+r[i][2]*r[i][2]) > radius_max);
    }
}

double rho_ws(double r)
{
    double r_ws = _R0_ * pow(_A_, 1./3);
    double rho0 = 3./4/M_PI/pow(_R0_,3);
    return rho0/(1+exp((r-r_ws)/_A_WS_));
}

void coords_generate(vector<vector<double> > &r, double radius_max)
{
    double rho0 = 3./4/M_PI/pow(_R0_,3);
    double radius, theta, phi;
    int nbr_points = r.size();

    /*for(int i=0 ; i<nbr_points ; i++)
    {
        do // Accepté avec proba P=rho(r)/rho(0)
        {
            radius = (radius_max+5*_A_WS_)*pow(alea(),1./3);
        }
        while(rho_ws(radius)/rho0 < alea());
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        r[i][0] = radius*sin(theta)*cos(phi);
        r[i][1] = radius*sin(theta)*sin(phi);
        r[i][2] = radius*cos(theta);
    }*/

    int counts = 0;
    for(int i=0 ; i<10*nbr_points ; i++)
    {
        radius = (radius_max+2*_A_WS_)*pow(alea(),1./3);
        if(rho_ws(radius)/rho0 < alea())
        {
            continue;
        }
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        r[counts][0] = radius*sin(theta)*cos(phi);
        r[counts][1] = radius*sin(theta)*sin(phi);
        r[counts][2] = radius*cos(theta);
        counts++;
        if(counts >= nbr_points)
        {
            break;
        }
    }
}

double fermi_momentum(double rho)
{
    return pow(3./2*M_PI*M_PI*rho,1./3) * _HBAR_C_;
}

void momenta_generate(vector<vector<double> > &r, vector<vector<double> > &p)
{
    double radius, theta, phi;
    int nbr_points = r.size();
    double pf;

    for(int i=0 ; i<nbr_points ; i++)
    {
        pf = fermi_momentum(rho_ws(module(r[i])));
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

double module(vector<double> &r)
{
    double sum = 0;
    for(int i=0 ; i<3 ; i++)
    {
        sum += r[i] * r[i];
    }
    return sqrt(sum);
}

double vect_module(vector<double> &r_real, vector<double> &r_box)
{
    double sum = 0;
    for(int i=0 ; i<3 ; i++)
    {
        sum += (r_real[i] - r_box[i]) * (r_real[i] - r_box[i]);
    }
    return sqrt(sum);
}

double gaussian(vector<double> &r_real, vector<double> &r_box)
{
    double r = vect_module(r_real, r_box);
    return 1/pow(sqrt(2*M_PI)*_SIGMA_,3)*exp(-r*r/2/_SIGMA_/_SIGMA_);
}

double simple_gaussian(double x1, double x2)
{
    double r = x1 - x2;
    return 1/sqrt(2*M_PI)/_SIGMA_*exp(-r*r/2/_SIGMA_/_SIGMA_);
}

double U(double rho)
{
    double rho0 = 3./4/M_PI/pow(_R0_,3); //fm-3
    return -356*rho/rho0 + 303*pow(rho/rho0,7./6);
}

int key2(int x, int y, int z)
{
    return x*_BOX_NBR_Y_*_BOX_NBR_Z_ + y*_BOX_NBR_Z_ + z;
}

void rho(vector<double> &rho_map, vector<vector<double> > &coords)
{
    //Initialize some variables
    double l0 = _L0_;
    double x1, y1, z1;
    int coords_size = coords.size();
    int rho_map_size = rho_map.size();
    int nbr_cells = floor(_SIGMA_NBR_*_SIGMA_/l0);
    vector<int> box_init(3);
    vector<vector<double> > r_box(2*nbr_cells+1, vector<double>(3));
    vector<vector<double> > gaus(2*nbr_cells+1, vector<double>(3));

    //Set rho map to zero
    for(int i=0 ; i<rho_map_size ; i++)
    {
        rho_map[i] = 0;
    }
    
    //Loop over test particles
    for(int i=0 ; i<coords_size ; i++)
    {
        //Set origin to particle coordinates
        for(int j=0 ; j<3 ; j++)
        {
            box_init[j] = floor(coords[i][j]/l0 + 0.5);
        }

        //Get position on cells and gaussian values
        for(int x=-nbr_cells ; x<=nbr_cells ; x++)
        {
            for(int j=0 ; j<3 ; j++)
            {
                r_box[x+nbr_cells][j] = (box_init[j] + x)*l0;
                gaus[x+nbr_cells][j] = simple_gaussian(r_box[x+nbr_cells][j], coords[i][j]);
            }
        }

        //Loop over considered cells
        for(int x=-nbr_cells ; x<=nbr_cells ; x++)
        {
            for(int y=-nbr_cells ; y<=nbr_cells ; y++)
            {
                for(int z=-nbr_cells ; z<=nbr_cells ; z++)
                {
                    //Coordinates on the grid
                    // _BOX_NBR_/2 to go on the middle of the grid
                    // +_BOX_NBR_ to avoid segmentation fault (always positive)
                    // %_BOX_NBR_ to avoid segmentation fault (always in the proper range)
                    x1 = (box_init[0] + x + _BOX_NBR_X_/2 + _BOX_NBR_X_) % _BOX_NBR_X_;
                    y1 = (box_init[1] + y + _BOX_NBR_Y_/2 + _BOX_NBR_Y_) % _BOX_NBR_Y_;
                    z1 = (box_init[2] + z + _BOX_NBR_Z_/2 + _BOX_NBR_Z_) % _BOX_NBR_Z_;

                    rho_map[key2(x1,y1,z1)] += gaus[x+nbr_cells][0]*gaus[y+nbr_cells][1]*gaus[z+nbr_cells][2];
                }
            }
        }
    }

    //Divide by _N_
    for(int i=0 ; i<rho_map_size ; i++)
    {
        rho_map[i] /= _N_;
    }
}

void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r)
{
    //Define some useful variables
    double l0 = _L0_;
    double gaus_tot;
    double pot;
    int box_init[3];
    int nbr_cells = floor(_SIGMA_NBR_*_SIGMA_/l0);
    int x1, y1, z1;
    vector<vector<double> > r_box(2*nbr_cells+1, vector<double>(3));
    vector<vector<double> > gaus(2*nbr_cells+1, vector<double>(3));

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        box_init[i] = floor(r[i]/l0 + 0.5);
        gradu[i] = 0;
    }

    //Get position on cells and gaussian values
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        for(int j=0 ; j<3 ; j++)
        {
            r_box[x+nbr_cells][j] = (box_init[j] + x)*l0;
            gaus[x+nbr_cells][j] = simple_gaussian(r[j], r_box[x+nbr_cells][j]);
        }
    }

    //Loop over considered cells
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        for(int y=-nbr_cells ; y<=nbr_cells ; y++)
        {
            for(int z=-nbr_cells ; z<=nbr_cells ; z++)
            {
                gaus_tot = gaus[x+nbr_cells][0]*gaus[y+nbr_cells][1]*gaus[z+nbr_cells][2];

                x1 = (box_init[0] + x + _BOX_NBR_X_/2 + _BOX_NBR_X_) % _BOX_NBR_X_;
                y1 = (box_init[1] + y + _BOX_NBR_Y_/2 + _BOX_NBR_Y_) % _BOX_NBR_Y_;
                z1 = (box_init[2] + z + _BOX_NBR_Z_/2 + _BOX_NBR_Z_) % _BOX_NBR_Z_;

                pot = U(rho_map[key2(x1,y1,z1)]);

                gradu[0] += (r[0] - r_box[x+nbr_cells][0]) * gaus_tot * pot;
                gradu[1] += (r[1] - r_box[y+nbr_cells][1]) * gaus_tot * pot;
                gradu[2] += (r[2] - r_box[z+nbr_cells][2]) * gaus_tot * pot;
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[i] *= l0*l0*l0/_SIGMA_/_SIGMA_;
    }
}

void minus_gradU2(vector<vector<double> > &gradu, vector<double> &rho_map, vector<vector<double> > &r)
{
    //Define some useful variables
    double l0 = _L0_;
    double gaus_tot;
    double pot;
    int box_init[3];
    int nbr_cells = floor(_SIGMA_NBR_*_SIGMA_/l0);
    int x1, y1, z1;
    vector<vector<double> > r_box(2*nbr_cells+1, vector<double>(3));
    vector<vector<double> > gaus(2*nbr_cells+1, vector<double>(3));

    //Loop over particles
    int part_nbr = r.size();
    for(int ipart=0 ; ipart<part_nbr ; ipart++)
    {

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        box_init[i] = floor(r[ipart][i]/l0 + 0.5);
        gradu[ipart][i] = 0;
    }

    //Get position on cells and gaussian values
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        for(int j=0 ; j<3 ; j++)
        {
            r_box[x+nbr_cells][j] = (box_init[j] + x)*l0;
            gaus[x+nbr_cells][j] = simple_gaussian(r[ipart][j], r_box[x+nbr_cells][j]);
        }
    }

    //Loop over considered cells
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        for(int y=-nbr_cells ; y<=nbr_cells ; y++)
        {
            for(int z=-nbr_cells ; z<=nbr_cells ; z++)
            {
                gaus_tot = gaus[x+nbr_cells][0]*gaus[y+nbr_cells][1]*gaus[z+nbr_cells][2];

                x1 = (box_init[0] + x + _BOX_NBR_X_/2 + _BOX_NBR_X_) % _BOX_NBR_X_;
                y1 = (box_init[1] + y + _BOX_NBR_Y_/2 + _BOX_NBR_Y_) % _BOX_NBR_Y_;
                z1 = (box_init[2] + z + _BOX_NBR_Z_/2 + _BOX_NBR_Z_) % _BOX_NBR_Z_;

                pot = U(rho_map[key2(x1,y1,z1)]);

                gradu[ipart][0] += (r[ipart][0] - r_box[x+nbr_cells][0]) * gaus_tot * pot;
                gradu[ipart][1] += (r[ipart][1] - r_box[y+nbr_cells][1]) * gaus_tot * pot;
                gradu[ipart][2] += (r[ipart][2] - r_box[z+nbr_cells][2]) * gaus_tot * pot;
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[ipart][i] *= l0*l0*l0/_SIGMA_/_SIGMA_;
    }

    }
}

double get_ubar(vector<double> &rho_map, vector<double> &r)
{
    //Define some useful variables
    double l0 = _L0_;
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    int nbr_cells = floor(_SIGMA_NBR_*_SIGMA_/l0);
    double ubar;
    int x1, y1, z1;

    //Initialize box coordinates
    ubar = 0;
    for(int i=0 ; i<3 ; i++)
    {
        box_init[i] = floor(r[i]/l0 + 0.5);
    }

    //Loop over considered cells
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        box[0] = (box_init[0] + x)*l0;
        for(int y=-nbr_cells ; y<=nbr_cells ; y++)
        {
            box[1] = (box_init[1] + y)*l0;
            for(int z=-nbr_cells ; z<=nbr_cells ; z++)
            {
                box[2] = (box_init[2] + z)*l0;
                gaus = gaussian(r, box);

                //(... + box_nbr)%box_nbr; to avoid negative number
                x1 = (box_init[0] + x + _BOX_NBR_/2 + _BOX_NBR_) % _BOX_NBR_;
                y1 = (box_init[1] + y + _BOX_NBR_/2 + _BOX_NBR_) % _BOX_NBR_;
                z1 = (box_init[2] + z + _BOX_NBR_/2 + _BOX_NBR_) % _BOX_NBR_;
                rho = rho_map[key(x1,y1,z1,_BOX_NBR_)];

                ubar += gaus * U(rho);
            }
        }
    }

    return ubar*l0*l0*l0;
}



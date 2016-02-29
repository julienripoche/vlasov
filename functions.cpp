#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

double alea()
{
    return (double)rand()/RAND_MAX;
}

vector<vector<double> > sphere(double radius_max, unsigned int nbr_points)
{
    double radius, theta, phi;
    vector<vector<double> > cart_coords(nbr_points, vector<double>(3));

    for(unsigned int i=0 ; i<nbr_points ; i++)
    {
        radius = radius_max*pow(alea(),1./3);
        theta = acos(1-2*alea());
        phi= 2*M_PI*alea();
        cart_coords[i][0] = radius*sin(theta)*cos(phi);
        cart_coords[i][1] = radius*sin(theta)*sin(phi);
        cart_coords[i][2] = radius*cos(theta);
    }

    return cart_coords;
}

void write(string const nomFichier, vector<vector<double> > &data)
{
    ofstream monFlux(nomFichier.c_str());
    unsigned int data_size = data.size();
    if(monFlux)
    {
        for(unsigned int i=0 ; i<data_size ; i++)
        {
            monFlux << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;
        }
    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier en écriture." << endl;
    }
}

void read(string const nomFichier, vector<vector<double> > &data)
{
    ifstream monFlux(nomFichier.c_str());
    vector<double> r(3);
    while(monFlux >> r[0] >> r[1] >> r[2])
    {
        data.push_back(r);
    }
}

/* Classical particle
double rho(double N, double L)
{
    double density = 0;
    double omega[3];
    omega[0] = 0.25;
    omega[1] = 0.5;
    omega[2] = 0.25;
    for(int i=-1 ; i<=1 ; i++)
    {
        for(int j=-1 ; j<=1 ; j++)
        {
            for(int k=-1 ; k<=1 ; k++)
            {
                density += omega[i+1]*omega[j+1]*omega[k+1]*;
            }
        }
    }
    return 1./N/L/L/L*
}
*/

double module(vector<double> &r_real, vector<double> &r_box)
{
    double sqr_sum = 0;
    for(unsigned int i=0 ; i<3 ; i++)
    {
        sqr_sum += (r_real[i] - r_box[i]) * (r_real[i] - r_box[i]);
    }
    return sqrt(sqr_sum);
}

double gaussian(vector<double> &r_real, vector<double> &r_box, double sigma)
{
    double r = module(r_real, r_box);
    return 1/pow(sqrt(2*M_PI)*sigma,3)*exp(-r*r/2/sigma/sigma);
}

unsigned int key(unsigned int x, unsigned int y, unsigned int z, unsigned int N)
{
    return x*N*N + y*N + z;
}

void rho_map(vector<double> &rho, vector<vector<double> > &coords, double sigma, double L, int box_size)
{
    int NA = coords.size();
    double l0 = double(L)/box_size;
    cout << "l0 = " << l0 << endl;
    vector<double> r(3);
    for(int x=0 ; x<box_size ; x++)
    {
        r[0] = (x-box_size/2.)*l0;
        for(int y=0 ; y<box_size ; y++)
        {
            r[1] = (y-box_size/2.)*l0;
            for(unsigned int z=0 ; z<box_size ; z++)
            {
                r[2] = (z-box_size/2.)*l0;
                for(unsigned int i=0 ; i<NA ; i++)
                {
                    rho[key(x,y,z,box_size)] += gaussian(r, coords[i], sigma);
                }
            }
        }
    }
}

double U(double rho)
{
    double r0 = 1.12; //fm
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    return -356*rho/rho0 + 303*pow(rho/rho0,7./6);
}

vector<double> minus_gradU(vector<double> &r, vector<double> &rho_map, double nbr_sigma, int N, double L, double sigma)
{
    //Define some useful variables
    double l0 = L/N;
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    vector<double> gradu(3,0);
    int nbr_cells = floor(nbr_sigma*sigma/l0);

    int x1, y1, z1;

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        //box_init[i] = (floor(r[i]/l0) + 0.5)*l0;
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
                gaus = gaussian(r, box, sigma);
                cout << gaus << " ";
                cout << key(x,y,z,N) << " ";

                x1 = (box_init[0] + x)%N;
                y1 = (box_init[1] + y)%N;
                z1 = (box_init[2] + z)%N;

                rho = rho_map[key(x1,y1,z1,N)];

                gradu[0] += (r[0] - (box_init[0]+x1+0.5)*l0) * gaus * U(rho);
                gradu[1] += (r[1] - (box_init[1]+y1+0.5)*l0) * gaus * U(rho);
                gradu[2] += (r[2] - (box_init[2]+z1+0.5)*l0) * gaus * U(rho);

                /*
                for(int i=0 ; i<3 ; i++)
                {
                    rho = rho_map[key(x,y,z,N)];
                    gradu[i] += (r[i] - (box_init[i]+[i]+0.5)*l0) * gaus * U(rho);
                }
                */
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[i] *= l0*l0*l0/sigma/sigma;
    }

    return gradu;
}

vector<double> minus_gradU2(vector<double> &r, vector<double> &rho_map, double nbr_sigma, int N, double L, double sigma)
{
    //Define some useful variables
    double l0 = L/N;
    double gaus;
    double rho;
    int box_init[3];
    vector<double> box(3);
    vector<double> gradu(3,0);
    int nbr_cells = floor(nbr_sigma*sigma/l0);

    int x1, y1, z1;

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        //box_init[i] = (floor(r[i]/l0) + 0.5)*l0;
        box_init[i] = floor(r[i]/l0);
    }

    return gradu;
}


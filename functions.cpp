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

double module(vector<double> coords)
{
    double sqr_sum = 0;
    for(int i=0 ; i<3 ; i++)
    {
        sqr_sum += coords[i]*coords[i];
    }
    return sqrt(sqr_sum);
}

double gaussian(vector<double> coords, double sigma)
{
    double modulus = module(coords);
    return 1/pow(sqrt(2*M_PI)*sigma,3)*exp(-modulus*modulus/2/sigma/sigma);
}

vector<double> gradU(vector<double> coords)
{
    //nbr_sigma sigma L U have to be arguments
    unsigned int nbr_sigma = 3;
    double L = 5;
    unsigned int N = 100;
    double l0 = L/N;
    double sigma = 0.05;
    double gaus;
    double U = -356.+303;
    double box_coords[3];
    vector<double> position_difference(3);
    vector<double> gradu(3,0);
    int nbr_cells = floor(nbr_sigma*sigma/l0);

    //Initialize box coordinates
    for(int i=0 ; i<3 ; i++)
    {
        box_coords[i] = floor(coords[i]/l0) + 0.5;
        cout << coords[i] << " " << box_coords[i] << endl;
    }

    //Loop over considered cells
    for(int x=-nbr_cells ; x<=nbr_cells ; x++)
    {
        position_difference[0] = coords[0] - (x + box_coords[0])*l0;
        for(int y=-nbr_cells ; y<=nbr_cells ; y++)
        {
            position_difference[1] = coords[1] - (y + box_coords[1])*l0;
            for(int z=-nbr_cells ; z<=nbr_cells ; z++)
            {
                position_difference[2] = coords[2] - (z + box_coords[2])*l0;

                gaus = gaussian(position_difference, sigma);

                for(int i=0 ; i<3 ; i++)
                {
                    gradu[i] += position_difference[i]*gaus*U;
                }
            }
        }
    }

    for(int i=0 ; i<3 ; i++)
    {
        gradu[i] *= -l0*l0*l0/sigma/sigma;
    }

    return gradu;
}



#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

double get_N(double a)
{
    double r0 = 1.12; //fm
    double R = r0 * pow(_A_, 1./3);
    double l = R/50;
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    double sum = 0;

    for(double r=0 ; r<2*R ; r+=l)
    {
        sum += r * r * rho0/(1+exp((r-R)/a));
    }

    sum *= 4*M_PI*l;

    return _NA_/sum;
}

int main()
{
    double r0 = 1.12;
    double R = r0 * pow(_A_, 1./3);
    double l = R/50;


    cout << "N new = " << get_N(0.5) << endl;

    ofstream wsFile("test_ws.gnu");
    for(double r=0 ; r<2*R ; r+=l)
    {
        wsFile << r << " " << rho_ws(r) << endl;
    }

    return 0;
}

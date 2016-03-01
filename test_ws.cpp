#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

double density(double r, double a, int A)
{
    double r0 = 1.12; //fm
    double R = r0 * pow(A, 1./3);
    double rho0 = 3./4/M_PI/pow(r0,3); //fm-3
    return rho0/(1+exp((r-R)/a));
}

using namespace std;

void f(double tab[][10])
{
    cout << tab[5][5] << endl;
}


int main()
{
    double a = 1.;
    int A = 56;
    double r0 = 1.12; //fm
    double R = r0 * pow(A, 1./3);
    double sum = 0;
    double l0 = R/50;

    ofstream wsFile("test_ws.gnu");
    for(double r=0 ; r<2*R ; r+=l0)
    {
        wsFile << r << " " << density(r,a,A) << endl;
        sum += r * r * density(r,a,A);
    }

    sum *= 4*M_PI*l0;

    cout << "sum = " << sum << endl;

    return 0;
}

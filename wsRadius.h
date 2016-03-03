#ifndef WSRADIUS
#define WSRADIUS

#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

class WsRadiusFinder
{
public:
    WsRadiusFinder();
    double rho_ws(double r);
    double rho_integration();
    double run();

private:
    double a_ws;
    double r_ws;
    double R;
    double rho0;
};

#endif

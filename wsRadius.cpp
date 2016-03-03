#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"
#include "wsRadius.h"

using namespace std;

WsRadiusFinder::WsRadiusFinder()
{
    this->a_ws = _A_WS_;
    this->r_ws = 0;
    this->R = _R0_ * pow(_A_,1./3);
    this->rho0 = 3./4/M_PI/pow(_R0_,3);
}

double WsRadiusFinder::rho_ws(double r)
{
    return this->rho0/(1+exp((r-this->r_ws)/this->a_ws));
}

double WsRadiusFinder::rho_integration()
{
    double step = this->a_ws/10;
    double r_max = this->r_ws + 5 * this->a_ws;
    double sum;

    for(double r=0 ; r<r_max ; r+=step)
    {
        sum += r * r * this->rho_ws(r);
    }

    sum *= 4*M_PI*step;

    return sum;
}

double WsRadiusFinder::run()
{
    double r_min = 0;
    double r_max = this->R;
    this->r_ws = (r_min + r_max)/2;
    double integral = this->rho_integration();

    cout << integral << " " << this->r_ws << endl;

    while(abs(integral-_A_)/_A_ > 1e-5)
    {
        if(integral > _A_)
        {
            r_max = this->r_ws;
            this->r_ws = (this->r_ws + r_min)/2;
            cout << endl;
        }
        else
        {
            r_min = this->r_ws;
            this->r_ws = (this->r_ws + r_max)/2;
            cout << endl;
        }
        integral = this->rho_integration();

        cout << integral << " " << this->r_ws << endl;
    }

    return this->r_ws;
}


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

int main()
{
    //double a_ws = _SIGMA_;
    WsRadiusFinder WsRadius;//(a_ws);

    //cout << WsRadius.rho_ws(6) << endl;
    //WsRadius.rho_integration();
    cout << "radius = " << WsRadius.run() << endl;

    return 0;
}

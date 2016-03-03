#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

int main()
{
    WsRadiusFinder WsRadius(0.001, 5);

    cout << WsRadius.rho_ws(0) << endl;

    //WsRadius.rho_integration();

    return 0;
}

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
    vector<double> r(3);
    vector<double> p(3);

    for(int i=0 ; i<3 ; i++)
    {
        r[i] = 1;
        p[i] = 2;
    }

    cout << "gaussian = " << gaussian(r, p, 1.) << endl;

    return 0;
}

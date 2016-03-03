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
    cout << key(0,0,0,100) << " " << key(0,0,15,100) << " " << key(0,15,0,100) << " " << key(15,0,0,100) << " " << key(15,15,15,100) << endl;

    return 0;
}

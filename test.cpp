#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include "functions.h"

using namespace std;

int main()
{
    srand48(time(NULL));

    cout << 1%25 << " " << 57%25 << " " << -2%25 << endl;

    cout << drand48() << endl;

    return 0;
}

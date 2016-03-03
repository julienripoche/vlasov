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
    srand(time(NULL));
    int N = 10000;
    double sum = 0;
    for(int i=0 ; i<N ; i++)
    {
        sum += alea();
        cout << alea() << endl;
    }
    cout << "mean = " << sum/N << endl;

    return 0;
}

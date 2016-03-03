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
    int i = 1;
    char fileName[50];
    sprintf(fileName, "file%d", i);
    cout << fileName << endl;

    string begin("aaa");
    string end("bbb");

    cout << begin + string(i) + end << endl;

    return 0;
}

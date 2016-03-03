#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std;

int main()
{
    int number;
    for(int i=0 ; i<10 ; i++)
    {
        cin >> number;
        cout << number%10 << endl;
    } 
    return 0;
}

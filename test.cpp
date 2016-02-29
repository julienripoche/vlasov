#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
//#include "functions.h"

using namespace std;

void dfunction(vector<vector<int> > &vi)
{
    vector<int> vect(2,3);
    vi.push_back(vect);
    vi.push_back(vect);
    vi.push_back(vect);
    vi.push_back(vect);
    vi.push_back(vect);
}


void ma_fonction(vector<int> &vi)
{
    vi.push_back(12);
}

int main()
{
    vector<int> vect(2,3);
    ma_fonction(vect);
    ma_fonction(vect);
    ma_fonction(vect);
    for(unsigned int i=0 ; i<vect.size() ; i++)
    {
        cout << vect[i] << endl;
    }

    cout << "2d" << endl;

    vector<vector<int> > vvv;
    vvv.push_back(vect);
    dfunction(vvv);
    for(unsigned int i=0 ; i<vvv.size() ; i++)
    {
        for(unsigned int j=0 ; j<vvv[i].size() ; j++)
        {
            cout << vvv[i][j] << " ";
        }
        cout << endl;
    }
}

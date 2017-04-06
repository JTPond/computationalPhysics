//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

//Function for returning cmath's random int 
int getRandomNumber(int min, int max){ 
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);  
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

int main(int argc, char const *argv[]) {
    //seeding the RNG with the current clock time
    srand(time(0));
    
    //constants for the LCG    
    int a = 63;
    int c = 319;
    int m = 65537;

    //seeds for the LCG
    int X0 = 561;
    int Y0 = 318;
    int Z0 = 654;

    //Both output files for rng's
    ofstream fs0,fs1;
    fs0.open("./dataFiles/bad_rng.dat");
    fs0.precision(16);
    fs1.open("./dataFiles/good_rng.dat");
    fs1.precision(16);

    //loop 
    for (int i = 0; i < 50000; i++) {
        //LCG calculating 
        X0 = (a*X0 + c)%m;
        Y0 = (a*Y0 + c)%m;
        Z0 = (a*Z0 + c)%m;

        //cmath RNG calcs
        int X1 = getRandomNumber(1,m-1);
        int Y1 = getRandomNumber(1,m-1);
        int Z1 = getRandomNumber(1,m-1);

        //outputs
        fs0 << X0 << "\t" << Y0 << "\t" << Z0 << "\n";
        fs1 << X1 << "\t" << Y1 << "\t" << Z1 << "\n";
    }
    return 0;
}

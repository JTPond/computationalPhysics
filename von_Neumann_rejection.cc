//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

int getRandomNumber(int min, int max){ 
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);  
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

double semi_circle(double x){
    double N = (2.0/(PI));
    if(abs(x) <= 1.0)return N*sqrt(1.0 - x*x);
    else return 0.0;
}    

double gaussian(double x){
    double N = 1.0/sqrt(2.0*PI);
    double U = -(x*x*.5);
    return N*exp(U);
}

int main(int argc, char const *argv[]) {
    srand(time(0));    
    ofstream fs0,fs1;
    fs0.open("./dataFiles/a.dat");
    fs0.precision(16);
    fs1.open("./dataFiles/b.dat");
    fs1.precision(16);

    for (int i = 0; i < 5000000; i++){
        double X0 = (((double) getRandomNumber(0,80000000))/10000000.)-4.0;
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);
        
        double Y0 = semi_circle(X0);
        double Y1 = gaussian(X0);   
        if(X1 < Y0)fs0 << X0 << "\n";
        if(X1 < Y1)fs1 << X0 << "\n";
    }
    fs0.close();
    fs1.close();
    return 0;
}

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

double semi_circle(double x, double a){
    double N = (2.0/(PI*a*a));
    if(abs(x) <= 1.0)return N*sqrt(a*a - x*x);
    else return 0.0;
}    

double brownian(double x){
    if(abs(x) < 3.0)return abs(x)/9.0;
    else return 0.0;
}

double get_random_start(){
    while(true){
        double X0 = (((double) getRandomNumber(0,60000000))/10000000.)-3.0;
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);    
        double Y0 = brownian(X0);
        if(X1 < Y0)return X0;
    }
}

double get_random_length(double a){
    while(true){
        double X0 = (((double) getRandomNumber(0,20000000))/10000000.)-1.0;
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);    
        double Y0 = semi_circle(X0,a);
        if(X1 < Y0)return X0;
    }
}

void RW_int(double a, long len){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/d.dat");
    double n_max = 1000000;
    fs.precision(8);
    //loop over possible values of N
    for(int n=1; n<=n_max; n++){ 
        double sum = get_random_start();
        for (long i = 0; i < len; i += 1){
            double x = get_random_length(a);
            sum += x;
        }
        fs << n << '\t' << sum << endl;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    RW_int(.5,500);
    return 0;
}



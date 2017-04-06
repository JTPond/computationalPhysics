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

double box(double x){
    if(abs(x) <= 1.0)return 0.5;
    else return 0.0;
}    

double dBump(double x){
    return 1.0*pow(x,4.0)*exp(-pow(x/0.75,2.0));
}

double get_random_length(double now){
    double X0 = (((double) getRandomNumber(0,24000000))/10000000.)-1.2;
    double X1 = (((double) getRandomNumber(0,10000000))/10000000.); 
    double Y0 = (dBump(now+X0)/dBump(now));
    if(X1 < Y0)return X0;
    else return 0.0;
    
}

void RW_int(double start){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    double n_max = 1000000;
    fs.precision(8);
    //loop over possible values of N
    for(int n=1; n<=n_max; n++){ 
        double sum = get_random_length(start);
        for (int i = 0; i < 100; i += 1){
            double x = get_random_length(sum);
            sum += x;
        }
        fs << n << '\t' << sum << endl;
    }
    fs.close();
}


int main(int argc, char const *argv[]) {
    srand(time(0));
    RW_int(0.0);
    return 0;
}



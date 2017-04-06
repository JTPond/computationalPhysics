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

double get_random_length(double a){
    while(true){
        double X0 = (((double) getRandomNumber(0,20000000))/10000000.)-1.0;
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);    
        double Y0 = semi_circle(X0,a);
        if(X1 < Y0)return X0;
    }
}

void RW_int(double a){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/g.dat");
    double n_max = 1000000;
    fs.precision(8);
    //loop over possible values of N
    for(int n=1; n<=n_max; n++){ 
        double sum = 0.0;
        for (int i = 0; i < 100; i += 1){
            double x = get_random_length(a);
            sum += x;
        }
        fs << n << '\t' << sum << endl;
    }
    fs.close();
}

double sigmax(double a){
    double sum11 = 0.0;
    for (int j = 0; j < 1000; j += 1){
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = 0; i < 100; i += 1){
            double x = get_random_length(a);
            sum1+=x;
            sum2+=x*x;
        }
        sum11+=sqrt((sum2/100.0) - (sum1*sum1/10000.0));
    }
    return sum11/100.0;
}

void f_sig(void){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/h.dat");
    fs.precision(8);
    for (double i = 0.1; i < .9; i += .2){
        fs << i << '\t' << sigmax(i) << endl;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    RW_int(.3);
    f_sig();
    return 0;
}



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

double EL(double x, double alpha){
    return alpha + x*x*(0.5 - 2.0*alpha*alpha);
}

double dBump(double x, double alpha){
    double norm = 1.0/sqrt(PI*(1.0/(alpha)));
    return pow(exp(-alpha*x*x),2.0);
}

double get_random_length(double now, double alpha){
    double X0 = (((double) getRandomNumber(0,20000000))/10000000.)-1.0;
    double X1 = (((double) getRandomNumber(0,10000000))/10000000.); 
    double Y0 = (dBump(now+X0,alpha)/dBump(now,alpha));
    if(X1 < Y0)return X0;
    else return 0.0;
    
}

void RW_int(double start){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    double n_max = 10000;
    fs.precision(8);
    //loop over possible values of alpha
    for (double alpha = 0.1; alpha < 1.8; alpha += 0.01){
        double Sum = 0.0;
        double Sum2 = 0.0;
        for(int n=1; n<=n_max; n++){ 
            double sum = get_random_length(start,alpha);
            for (int i = 0; i < 1000; i += 1){
                double x = get_random_length(sum,alpha);
                sum += x;
            }
            Sum += EL(sum,alpha);
            Sum2 += pow(EL(sum,alpha),2.0);
        }
        fs << alpha << "\t" << Sum/n_max << "\t" << Sum2/n_max - pow(Sum/n_max,2.0)<<"\n";
    }
    fs.close();
}

void descent(double start){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    double n_max = 10000;
    double gamma = 0.05;
    double alpha = start;
    fs.precision(8);
    for (int N = 0; N < 100; N += 1){
        fs << N << "\t" << alpha <<"\n";
        double Sum = 0.0;
        double Sum2 = 0.0;
        double Sum3 = 0.0;
        for(int n=1; n<=n_max; n++){ 
            double sum = get_random_length(0.0,alpha);
            for (int i = 0; i < 1000; i += 1){
                double x = get_random_length(sum,alpha);
                sum += x;
            }
            Sum -= EL(sum,alpha)*sum*sum;
            Sum2 += sum*sum;
            Sum3 += sum;
        }
        alpha -= gamma*2.0*(Sum/n_max + EL(Sum3/n_max,alpha)*Sum2/n_max);
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    //RW_int(0.0);
    descent(0.7);
    return 0;
}



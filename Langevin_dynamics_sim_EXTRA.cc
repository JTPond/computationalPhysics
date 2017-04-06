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

double force(double gamma, double dt, double KbT, double x){
    double sigma2 = 2.0*gamma*KbT*dt;
    double norm = 1.0/(sqrt(2.0*sigma2*PI));
    return exp((-(x*x))/(2*sigma2));
}

double random_force(double gamma, double dt, double KbT){
    while(true){
        double X0 = (((double) getRandomNumber(0,60000000))/10000000.) - 3.0;
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);
        double Y0 = force(gamma,dt,KbT,X0);
        if (X1 < Y0){return X0;}
    }
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    double dt = 0.02;
    double gamma = 0.1;
    double KbT = 0.5;
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/f.dat");
    double x[10000]; double v[10000];
    for (size_t i = 0; i < 10000; i++) {
        x[i] = 0.0; v[i] = 0.0;
    }
    for (size_t t = 0; t < 3500/dt; t++) {
        for (size_t j = 0; j < 10000; j++) {
            double alpha = exp(-1.0*gamma*dt);
            double fn = random_force(gamma,dt,KbT);
            x[j] = x[j] + v[j]*dt;
            v[j] = alpha*v[j] + (fn);
        }
        double X = 0.0;
        for (size_t k = 0; k < 10000; k++) {
            X+=x[k]*x[k];
        }
        fs << t*dt << '\t' << X/10000.0 << '\n';
    }
    fs.close();
    return 0;
}

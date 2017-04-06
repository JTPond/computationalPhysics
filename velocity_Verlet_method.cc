//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double potential(double r){
    return pow(1.0/r,12.0) - 2*pow(1.0/r,6.0);
}

double force(double r){
    return 12.0*(pow(1.0/r,13.0) - pow(1.0/r,7));
}

int main(int argc, char const *argv[]) {
    double dt = 0.008;
    double x1 = -2.0; double x2 = 2.0;
    double v1 =  0.0; double v2 = 0.0;
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    for (size_t i = 0; i < 250.0/dt; i++) {
        double K = (v1*v1)/2.0 + (v2*v2)/2.0;
        double P = potential(abs(x2-x1));
        fs << i*dt << '\t' << x1 << '\t' << x2 << '\t' << v1 << '\t' << v2 << '\t' << K << '\t' << P << '\n';
/*For epsilon. Note: Lines 28 and 30 are never live at the same time*/
        fs << i*dt << '\t' << (K+P) + 0.000488221645891 << '\n';

        double fn = force(abs(x2-x1));
        //Notice the Negative Force term on x1 and v1
        x1 = x1 + v1*dt - 0.5*fn*dt*dt; x2 = x2 + v2*dt + 0.5*fn*dt*dt;
        v1 = v1 - 0.5*(fn + force(abs(x2-x1)))*dt; v2 = v2 + 0.5*(fn + force(abs(x2-x1)))*dt;

    }
    fs.close();

    return 0;
}

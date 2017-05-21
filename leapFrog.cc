//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

double dt =  0.0000002;

const int N = 200;

double big_U[N];
double big_U2[N];

void set_string(){
    for (int i = 0; i < N; i+=1.0) {
        double x = i/100.0;
        double xm1 = (i-1)/100.0;
        double xm2 = (i-2)/100.0;
        double xp1 = (i+1)/100.0;
        double xp2 = (i+2)/100.0;
        big_U[i] = cos(PI*x) - ((dt)/(0.01))*(cos(PI*xp1) + cos(PI*x) + cos(PI*xm1))*(cos(PI*xp1) - cos(PI*xm1)) - ((dt)/(2.0*pow(0.01,3.0)))*(cos(PI*xp2) - 2.0*cos(PI*xp1) + 2.0*cos(PI*xm1) - cos(PI*xm2));
        big_U2[i] = cos(PI*x);
    }
}

void reset_string() {
    double temp_U[N];

    for (int j = 0; j < N; j++)  temp_U[j] = big_U[j];

    big_U[0] = big_U2[0] - ((dt)/(3.0*0.01))*(temp_U[1] + temp_U[0] + temp_U[199])*(temp_U[1] - temp_U[199]) - ((pow(0.022,2.0)*dt)/pow(0.01,3.0))*(temp_U[2] - 2.0*temp_U[1] + 2.0*temp_U[199] - temp_U[198]);

    big_U[1] = big_U2[1] - ((dt)/(3.0*0.01))*(temp_U[2] + temp_U[1] + temp_U[0])*(temp_U[2] - temp_U[0]) - ((pow(0.022,2.0)*dt)/pow(0.01,3.0))*(temp_U[3] - 2.0*temp_U[2] + 2.0*temp_U[0] - temp_U[199]);

    for (int i = 2; i < N-2; i++) {
        big_U[i] = big_U2[i] - ((dt)/(3.0*0.01))*(temp_U[i+1] + temp_U[i] + temp_U[i-1])*(temp_U[i+1] - temp_U[i-1]) - ((pow(0.022,2.0)*dt)/pow(0.01,3.0))*(temp_U[i+2] - 2.0*temp_U[i+1] + 2.0*temp_U[i-1] - temp_U[i-2]);
    }

    big_U[198] = big_U2[198] - ((dt)/(3.0*0.01))*(temp_U[199] + temp_U[198] + temp_U[197])*(temp_U[199] - temp_U[197]) - ((pow(0.022,2.0)*dt)/pow(0.01,3.0))*(temp_U[0] - 2.0*temp_U[199] + 2.0*temp_U[197] - temp_U[196]);

    big_U[199] = big_U2[199] - ((dt)/(3.0*0.01))*(temp_U[0] + temp_U[199] + temp_U[198])*(temp_U[0] - temp_U[198]) - ((pow(0.022,2.0)*dt)/pow(0.01,3.0))*(temp_U[1] - 2.0*temp_U[0] + 2.0*temp_U[198] - temp_U[197]);

    for (int k = 0; k < N; k++)  big_U2[k] = temp_U[k];
}

void run_leapFrog(){
    set_string();
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    double count = 0.0;
    for (double i = 0; i < 3.5/dt; i++) {
        if (count == 0.0){
            cout << i*dt << "\n";
            fs << i*dt;
            for (size_t j = 0; j < N; j+=1.0) {
                fs << ";" << big_U[j];
            }
            fs << "\n";
        }
        reset_string();
        if (count < 50000) count +=1.0;
        else count = 0.0;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    run_leapFrog();
    return 0;
}

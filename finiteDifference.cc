//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

const int N = 700;

double big_U[N];

double energy(double dx){
    double kin = 0.0;
    for (size_t i = 0; i < N; i++) {
        kin += pow(big_U[i],2.0)*dx;
    }
    return kin;
}

void set_string(double w, double dx){
    for (int i = 0; i < N; i+=1.0) {
        double x = (i*dx) - 20.0;
        big_U[i] = exp(-pow(x,2.0)/pow(w,2.0));
    }
}

void reset_string(double dt, double dx, double v) {
    double temp_U[N];
    for (int j = 0; j < N; j++)  temp_U[j] = big_U[j];
    for (int i = 1; i < N-1; i++) {
        big_U[i] = 0.5*(temp_U[i+1] + temp_U[i-1]) - ((v*dt)/(2.0*dx))*(temp_U[i+1] - temp_U[i-1]);
    }
}

void run_energy_cons(){
    double dt = 0.14;
    double dx = 0.1;
    double w = 5.0;
    double v = 1.0;
    set_string(w,dx);
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/c.dat");
    for (double i = 0; i < (20.0)/dt; i++) {
        if(fmod(i*dt,5.0) < 0.1){
            fs << i*dt;
            for (size_t j = 0; j < N; j+=1.0) {
                double x = (j*dx) - 20.0;
                fs << ";" << x << "," << big_U[j];
            }
            fs << "\n";
        }
        reset_string(dt,dx,v);
    }
    fs.close();
}

void run_energy_cons_line(){
    double dx = 0.1;
    double w = 5.0;
    double v = 1.0;
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    for (double dt = 0.02; dt < 0.18; dt+=0.02){
        set_string(w,dx);
        double E0 = energy(dx);
        fs << dt;
        for (double i = 0; i < (20.0)/dt; i++) {
            fs << ";" << i*dt;
            fs << "," << energy(dx)/E0;
            reset_string(dt,dx,v);
        }
        fs << "\n";
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    run_energy_cons();
    run_energy_cons_line();
    return 0;
}

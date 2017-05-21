//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

const int N = 36;

double big_U[N+2];
double big_V[N+2];

double potential(double k){
    double pot = 0.0;
    for (size_t i = 0; i < N+1; i++) {
        pot += (k/2.0)*pow(big_U[i+1] - big_U[i],2.0);
    }
    return pot;
}

double kinetic(double m){
    double kin = 0.0;
    for (size_t i = 1; i < N+1; i++) {
        kin += (m/2.0)*pow(big_V[i],2.0);
    }
    return kin;
}

double E_local(int n, double m, double k){
    double pot = (k/4.0)*(pow(big_U[n+1] - big_U[n],2.0) + pow(big_U[n] - big_U[n-1],2.0));
    double kin = (m/2.0)*pow(big_V[n],2.0);
    return kin + pot;     
}

double force_n(int n, double k, double alpha){
    double kterm = k*(big_U[n+1] + big_U[n-1] - 2.0*big_U[n]);
    double alphaterm = alpha*(pow(big_U[n+1] - big_U[n],2.0) - pow(big_U[n] - big_U[n-1],2.0));
    return kterm + alphaterm;
}

double epsilon_n(double n, double m){
    double norm = sqrt(2.0/(N + 1));
    return norm*sin((m*n*PI)/(N+1));
}

double omega_m(double m, double k){
    double norm = 2.0*sqrt(k/m);
    return norm*sin((m*PI)/(2.0*(N+1)));
}

double E_m(double m, double k){
    double Qm = 0.0;
    double Q_dot_m = 0.0;
    for (size_t i = 1; i < N+1; i++) {
        double eps = epsilon_n((double) i,m);
        Qm += big_U[i]*eps;
        Q_dot_m += big_V[i]*eps;
    }
    return (m/2.0)*(pow(Q_dot_m,2.0) + pow(omega_m(m,k),2.0)*pow(Qm,2.0));
}

void set_string(double A, double B, double m, double sigma){
    big_U[0] = 0.0; big_U[N+1] = 0.0;
    big_V[0] = 0.0; big_V[N+1] = 0.0;
    for (size_t i = 1; i < N+1; i+=1.0) {
        big_U[i] = A*epsilon_n(i,m);
        big_V[i] = B*exp((-pow(i-((N+2)/2.0),2.0))/(sigma*sigma));
    }
}


void reset_string(double dt, double k, double alpha) {
    for (size_t i = 1; i < N+1; i++) {
        double fn1 = force_n(i, k, alpha);
        big_U[i] = big_U[i] + big_V[i]*dt + 0.5*fn1*dt*dt;
        double fn2 = force_n(i, k, alpha);
        big_V[i] = big_V[i] + 0.5*(fn1 + fn2)*dt;
    }
}

void run_energy_cons(){
    double dt = 0.07;
    double k = 1.0;
    double m = 1.0;
    double alpha = 0.0;
    double A = 0.2;
    double B = 0.0;
    double sigma = 1.0;
    set_string(A,B,m,sigma);
    double w = omega_m(m,k);
    double E0 = kinetic(m) + potential(k);
    double end = 2.0*2.0*PI/(w);
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/g.dat");
    for (size_t i = 0; i < end; i+=end/50.0) {
        for (size_t d = 0.0; d < i/dt; d++) reset_string(dt,k,alpha);
        double K = kinetic(m);
        double P = potential(k);
        fs << i*w/(2.0*PI) << '\t' << K << '\t' << P << '\t' << abs(E0 - (K+P))/E0 <<'\n';
        set_string(A,B,m,sigma);
    }
    fs.close();    
}

void run_energy_local(){
    double dt = 0.04;
    double k = 1.0;
    double m = 1.0;
    double alpha = -0.99;
    double A = 0.0;
    double B = 0.8;
    double sigma = 10.0;
    set_string(A,B,m,sigma);
    double w = omega_m(m,k);
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/n.dat");
    double end = 0.125*2.0*PI/(w);
    for (size_t i = 0.0; i < end; i+= end/5.0) {
        for (size_t d = 0.0; d < i/dt; d++) reset_string(dt,k,alpha);
        fs << i*w/(2.0*PI) << '\t'; 
        for (size_t j = 1; j < N+1; j++){
            fs << E_local(j,m,k) << '\t';
        }
        fs << "\n";
        set_string(A,B,m,sigma);

    }
    fs.close();    
}

void run_energy_mode(){
    double dt = 0.02;
    double k = 1.0;
    double m = 1.0;
    double alpha = 0.0;
    double A = 0.5;
    double B = 0.0;
    double sigma = 1.0;
    set_string(A,B,m,sigma);
    double w = omega_m(m,k);
    double end = 140.0*2.0*PI/(w);
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/m_alph0.dat");
    for (size_t i = 0.0; i < end; i+=end/50.0) {
        fs << i*w/(2.0*PI);
        for (double j = 1.0; j < 5.0; j += 1.0) fs << '\t' << E_m(j,k);
        fs << "\n";
        for (double j = 0.0; j < end/(25.0*dt); j += 1.0) reset_string(dt,k,alpha);
    }
    fs.close();    
}


int main(int argc, char const *argv[]) {
    run_energy_cons();
    run_energy_local();
    run_energy_mode();

    return 0;
}

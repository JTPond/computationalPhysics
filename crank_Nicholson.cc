//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;
double PI = 2. * asin(1.0);

double D = 2.0;
double w = 10.0;
double dx = 0.1;
double dt = 8.0;
const int xLen = 120.0/dx - 1;

double r = D*dt/(2.0*dx*dx);

mat big_co = mat(xLen,xLen);

mat big_co_l = mat(xLen,xLen);

vec big_U = vec(xLen);

void set_string(){
    big_co(0,0) = (1.0 + 2.0*r);
    big_co(1,0) = -r;
    for (int i = 1; i < xLen-1; i++){
        big_co(i+1,i) = -r;
        big_co(i,i) = (1.0 + 2.0*r);
        big_co(i-1,i) = -r;
    }
    big_co(xLen-1,xLen-1) = (1.0 + 2.0*r);
    big_co(xLen - 2,xLen-1) = -r;

    big_co_l(0,0) = (1.0 - 2.0*r);
    big_co_l(1,0) = r;
    for (int i = 1; i < xLen-1; i++){
        big_co_l(i+1,i) = r;
        big_co_l(i,i) = (1.0 - 2.0*r);
        big_co_l(i-1,i) = r;
    }
    big_co_l(xLen-1,xLen-1) = (1.0 - 2.0*r);
    big_co_l(xLen - 2,xLen-1) = r;

    for (int i = 0; i < xLen; i++){
        double x = (((double) i+1)*dx - (60.0));
        big_U(i) = exp(-(x*x)/(w*w));
    }
}
void reset_string(){
    //cout << big_co_l.n_rows << " , " << big_U.n_elem << "\n";
    vec big_d = big_co_l*big_U;

    big_co(1,0) = big_co(1,0)/big_co(0,0);
    for (int i = 1; i < xLen-1; i++){
        big_co(i+1,i) = big_co(i+1,i)/(big_co(i,i) - big_co(i-1,i)*big_co(i,i-1));
    }
    big_d(0) = big_d(0)/big_co(0,0);
    for (int i = 1; i < xLen; i++){
        big_d(i) = (big_d(i) - big_co(i-1,i)*big_d(i-1))/(big_co(i,i) - big_co(i-1,i)*big_co(i,i-1));
    }
    big_U(xLen-1) = big_d(xLen-1);
    for (int i = xLen-2; i >= 0; i--){
        big_U(i) = big_d(i) - big_co(i+1,i)*big_U(i+1);
    }

//    big_co_l = big_co;

    big_co(0,0) = (1.0 + 2.0*r);
    big_co(1,0) = -r;
    for (int i = 1; i < xLen-1; i++){
        big_co(i+1,i) = -r;
        big_co(i,i) = (1.0 + 2.0*r);
        big_co(i-1,i) = -r;
    }
    big_co(xLen-1,xLen-1) = (1.0 - 2.0*r);
    big_co(xLen - 2,xLen-1) = r;

    big_co_l(0,0) = (1.0 - 2.0*r);
    big_co_l(1,0) = r;
    for (int i = 1; i < xLen-1; i++){
        big_co_l(i+1,i) = r;
        big_co_l(i,i) = (1.0 - 2.0*r);
        big_co_l(i-1,i) = r;
    }
    big_co_l(xLen-1,xLen-1) = (1.0 - 2.0*r);
    big_co_l(xLen - 2,xLen-1) = r;
}

double analytical(double x,double t){
    double norm = 1.0/sqrt(1.0 + ((4.0*D*t)/(w*w)));
    return  norm * exp(-(x*x)/(w*w + 4.0*D*t));
}

void run_crankN_diff(){
    set_string();
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    double count = 0.0;
    for (double i = 0; i < 69.0/dt; i++) {
        if (count == 0.0){
            cout << i*dt << "\n";
            fs << i*dt;
            for (size_t j = 0; j < xLen; j+=1.0) {
                fs << ";" << big_U(j);
            }
            fs << "\n";
        }
        reset_string();
        if (count < 10/dt - 1.0) count +=1.0;
        else count = 0.0;
    }
    fs.close();
}

void run_crankN_error(){
    set_string();
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/f.dat");
    for (double i = 0; i < 200.0/dt; i++) {
        double a = analytical(0.0,i*dt);
        fs << i*dt << "\t"<< abs(a - big_U(599))/a <<"\n";
        reset_string();
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    run_crankN_error();
    return 0;
}

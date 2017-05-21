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

double ET = 3.0/2.0;

int M0 = 1000;

mat walkers = randn<mat>(3,1000);

double V(vec r){
    return accu(r%r)/2.0;
}

double W(vec rp){
    return exp(ET - V(rp));
}

double G(vec rp, vec r){
    double norm = 1.0/sqrt(2.0*PI);
    return norm*exp(-accu((rp - r)%(rp - r))/2.0);
}

void update_string(){
    int additions = 0;
    for (size_t i = 0; i < M0; i++) {
        walkers.col(i) += randn<vec>(3);
        double s;
        double w = W(walkers.col(i));
        double wfract = modf(w, &s);
        vec test = randu<vec>(1);
        if (test(0) < (w - s)) s += 1.0;
        if (s == 0){
             walkers.shed_col(i);
             i--;
             M0--;
        }
        else{
            walkers.resize(3,M0+additions+ ((int) s-1.0));
            for (size_t j = 0; j < ((int) s-1.0); j++) {
                walkers.col(M0+additions+ j) = walkers.col(i);
            }
            additions += ((int) s-1.0);
        }
    }
    M0 += additions;
    ET += 0.5*log(2000.0/((double) M0));
}

void run() {
    ofstream fs,fs1;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    fs1.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    double count = 0.0;
    while (count < 300.0) {
        update_string();
        cout << M0 << "\n";
        fs1 << M0 << "\t"<< ET <<"\n";
        count += 1.0;
    }
    walkers.t().print(fs);
    fs.close();
    fs1.close();
}

int main(int argc, char const *argv[]) {
    run();
    return 0;
}

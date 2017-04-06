//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

int big_X[1000][50];

int big_Y[1000][50];

double boltz(double E,double T){
    double k = 1.0;//8.617*pow(10,-5);
    return exp(-E/(k*T));
}

int getRandomNumber(int min, int max){
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

void set_string(){
    for (size_t j = 0; j < 1000; j++) {
        for (size_t i = 0; i < 50; i++) {
            int v = getRandomNumber(0,1);
            if(v == 0)v=-1;
            big_X[j][i] = v;
            big_Y[j][i] = v;
        }
    }
}

double calc_X_E(int k){
    double buff = 0.0;
    for (size_t i = 0; i < 50; i++) {
        if(i < 49){
            buff = buff + big_X[k][i]*big_X[k][i+1];
        }
        else{
            buff = buff + big_X[k][49]*big_X[k][0];
        }
    }
    return -buff;
}

double calc_Y_E(int k){
    double buff = 0.0;
    for (size_t i = 0; i < 50; i++) {
        if(i < 49){
            buff = buff + big_Y[k][i]*big_Y[k][i+1];
        }
        else{
            buff = buff + big_Y[k][49]*big_Y[k][0];
        }
    }
    return -buff;
}

double calc_avg_E(){
    double buff = 0.0;
    for (size_t k = 0; k < 1000; k++) {
        buff = buff + calc_X_E(k);
    }
    return buff/1000.0;
}

double calc_avg_E2(){
    double buff = 0.0;
    for (size_t k = 0; k < 1000; k++) {
        double cE = calc_X_E(k);
        buff = buff + cE*cE;
    }
    return buff/1000.0;
}

double calc_M(int k){
    double buff = 0.0;
    for (size_t i = 0; i < 50; i++) {
        buff = buff + big_X[k][i];
    }
    return -buff;
}

double calc_avg_M(){
    double buff = 0.0;
    for (size_t k = 0; k < 1000; k++) {
        buff = buff + calc_M(k);
    }
    return buff/1000.0;
}

double calc_avg_M2(){
    double buff = 0.0;
    for (size_t k = 0; k < 1000; k++) {
        double cM = calc_M(k);
        buff = buff + cM*cM;
    }
    return buff/1000.0;
}

double calc_epsilon(){
    return calc_avg_E()/50.;
}

double calc_c(double T){
    return (calc_avg_E2() - pow(calc_avg_E(),2))/(50.0 * T*T);
}

double calc_m(){
    return calc_avg_M()/50.;
}

double calc_chi(double T){
    return (calc_avg_M2() - pow(calc_avg_M(),2))/(50.0 * T);
}

void reset_string(double T){
    for (size_t i = 0; i < 1000; i++) {
        int j = getRandomNumber(0,49);
        big_Y[i][j] = big_Y[i][j] * -1;
        double del_E = calc_Y_E(i) - calc_X_E(i);
        double X0 = (((double) getRandomNumber(0,10000000))/10000000.);
        double B0 = boltz(del_E,T);
        if(del_E < 0){
            big_X[i][j] = big_X[i][j] * -1;
        }
        else if(X0 < B0){
            big_X[i][j] = big_X[i][j] * -1;
        }
        else{
            big_Y[i][j] = big_Y[i][j] * -1;
        }
    }
}

void RW_int(){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    for (double i = 0.1; i < 5.1; i += .1){
        for (size_t j = 0; j < 1000000; j++) {
            reset_string(i);
        }
        fs << i << "\t" <<  calc_epsilon() << "\t" << calc_c(i) << "\t" << calc_m() << "\t" << calc_chi(i) << endl;
        set_string();
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    set_string();
    RW_int();
    return 0;
}

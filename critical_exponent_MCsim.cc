//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

const int N = 50;

int big_X[1000][N][N];

double rateT = 0.0;
double rateB = 0.0;

double num_clusters = 0.0;
double clusterss = 0.0;

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
        for (size_t i = 0; i < N; i++) {
            for (size_t k = 0; k < N; k++) {
                int v = getRandomNumber(0,1);
                if(v == 0)v=-1;
                big_X[j][i][k] = v;
            }
        }
    }
}

double get_H_site(int i, int j, int k){
    int neighbors[2][2];
    if (j-1 < 0) neighbors[0][0] = (N-1);
    else neighbors[0][0] = j-1;
    if (j+1 > (N-1)) neighbors[0][1] = 0;
    else neighbors[0][1] = j+1;
    if (k-1 < 0) neighbors[1][0] = (N-1);
    else neighbors[1][0] = k-1;
    if (k+1 > (N-1)) neighbors[1][1] = 0;
    else neighbors[1][1] = k+1;
    double H = 0.0;
    H += (big_X[i][neighbors[0][0]][k]);
    H += (big_X[i][neighbors[0][1]][k]);
    H += (big_X[i][j][neighbors[1][0]]);
    H += (big_X[i][j][neighbors[1][1]]);
    return big_X[i][j][k]*H;
}

double get_H_hat(){
    double E = 0.0;
    for (size_t i = 0; i < 1000; i++) {
        double H = 0.0;
        bool gate = true;
        for (size_t k = 0; k < N; k++) {
            if (gate){
                for (size_t j = 0; j < N; j+=2) {
                    H += get_H_site(i,j,k);
                }
                gate = false;
            }
            else{
                for (size_t j = 1; j < N; j+=2) {
                    H += get_H_site(i,j,k);
                }
                gate = true;
            }
        }
        E -= H;
    }
    return E/1000.0;
}

double get_H_2_hat(){
    double E = 0.0;
    for (size_t i = 0; i < 1000; i++) {
        double H = 0.0;
        bool gate = true;
        for (size_t k = 0; k < N; k++) {
            if (gate){
                for (size_t j = 0; j < N; j+=2) {
                    H += get_H_site(i,j,k);
                }
                gate = false;
            }
            else{
                for (size_t j = 1; j < N; j+=2) {
                    H += get_H_site(i,j,k);
                }
                gate = true;
            }
        }
        E += H*H;
    }
    return E/(1000.0);
}

double get_C_hat(double T, double H_hat){
    return (get_H_2_hat() - pow(H_hat,2))/(T*T*N*N);
}

double get_M_pow(int p){
    double M;
    double m = 0.0;
    for (size_t i = 0; i < 1000; i++) {
        M = 0.0;
        for (size_t j = 0; j < N; j++) {
            for (size_t k = 0; k < N; k++) {
                M += big_X[i][j][k];
            }
        }
        m += abs(pow(M,p));
    }
    return m/(1000.0);
}

double get_m(){
    return get_M_pow(1)/(N*N);
}

double get_chi(double T){
    return (get_M_pow(2) - pow(get_M_pow(1),2.0))/(T*N*N);
}

double get_B4(){
    return 1.0 - (get_M_pow(4))/(3.0*pow(get_M_pow(2),2.0));
}

double get_delE(int i, int j, int k){
    int neighbors[2][2];
    if (j-1 < 0) neighbors[0][0] = (N-1);
    else neighbors[0][0] = j-1;
    if (j+1 > (N-1)) neighbors[0][1] = 0;
    else neighbors[0][1] = j+1;
    if (k-1 < 0) neighbors[1][0] = (N-1);
    else neighbors[1][0] = k-1;
    if (k+1 > (N-1)) neighbors[1][1] = 0;
    else neighbors[1][1] = k+1;
    double Hmol = 0.0;
    Hmol += (big_X[i][neighbors[0][0]][k]);
    Hmol += (big_X[i][neighbors[0][1]][k]);
    Hmol += (big_X[i][j][neighbors[1][0]]);
    Hmol += (big_X[i][j][neighbors[1][1]]);
    return 2.0*big_X[i][j][k]*Hmol;
}


bool in(int i, int j, double cluster_length, int cluster[N*N][2]){
    for (size_t k = 0; k < cluster_length; k++) {
        if (i == cluster[k][0] && j == cluster[k][1]) {
            return true;
        }
    }
    return false;
}

void reset_string_cluster(double T){
    for (size_t i = 0; i < 1000; i++) {
        int cluster_length = 1;
        int cluster[N*N][2];
        int f_length = 1;
        int f[N*N][2];
        int j = getRandomNumber(0,(N-1));
        int k =  getRandomNumber(0,(N-1));
        cluster[0][0] = j; cluster[0][1] = k;
        f[0][0] = j; f[0][1] = k;
        while (f_length != 0){
            int f_new_length = 0;
            int f_new[N*N][2];
            for (size_t t = 0; t < f_length; t++) {
                int i1 = f[t][0]; int i2 = f[t][1];

                int neighbors[2][2];
                if (i1-1 < 0) neighbors[0][0] = (N-1);
                else neighbors[0][0] = i1-1;
                if (i1+1 > (N-1)) neighbors[0][1] = 0;
                else neighbors[0][1] = i1+1;
                if (i2-1 < 0) neighbors[1][0] = (N-1);
                else neighbors[1][0] = i2-1;
                if (i2+1 > (N-1)) neighbors[1][1] = 0;
                else neighbors[1][1] = i2+1;

                if (big_X[i][neighbors[0][0]][i2] == big_X[i][i1][i2] && !in(neighbors[0][0],i2,cluster_length,cluster)){
                    double X0 = (rand() % 10000000)/10000000.0;
                    double B0 = 1.0-exp(-2.0/T);
                    if (X0 < B0) {
                        f_new[f_new_length][0] = neighbors[0][0]; f_new[f_new_length][1] = i2;
                        f_new_length++;
                        cluster[cluster_length][0] = neighbors[0][0]; cluster[cluster_length][1] = i2;
                        cluster_length++;
                    }
                }
                if (big_X[i][neighbors[0][1]][i2] == big_X[i][i1][i2] && !in(neighbors[0][1],i2,cluster_length,cluster)){
                    double X0 = (rand() % 10000000)/10000000.0;
                    double B0 = 1.0-exp(-2.0/T);
                    if (X0 < B0) {
                        f_new[f_new_length][0] = neighbors[0][1]; f_new[f_new_length][1] = i2;
                        f_new_length++;
                        cluster[cluster_length][0] = neighbors[0][1]; cluster[cluster_length][1] = i2;
                        cluster_length++;
                    }
                }
                if (big_X[i][i1][neighbors[1][0]] == big_X[i][i1][i2] && !in(i1,neighbors[1][0],cluster_length,cluster)){
                    double X0 = (rand() % 10000000)/10000000.0;
                    double B0 = 1.0-exp(-2.0/T);
                    if (X0 < B0) {
                        f_new[f_new_length][0] = i1; f_new[f_new_length][1] = neighbors[1][0];
                        f_new_length++;
                        cluster[cluster_length][0] = i1; cluster[cluster_length][1] = neighbors[1][0];
                        cluster_length++;
                    }
                }
                if (big_X[i][i1][neighbors[1][1]] == big_X[i][i1][i2] && !in(i1,neighbors[1][1],cluster_length,cluster)){
                    double X0 = (rand() % 10000000)/10000000.0;
                    double B0 = 1.0-exp(-2.0/T);
                    if (X0 < B0) {
                        f_new[f_new_length][0] = i1; f_new[f_new_length][1] = neighbors[1][1];
                        f_new_length++;
                        cluster[cluster_length][0] = i1; cluster[cluster_length][1] = neighbors[1][1];
                        cluster_length++;
                    }
                }
            }
            f_length = f_new_length;
            for (size_t z = 0; z < f_new_length; z++) {
                f[z][0] = f_new[z][0]; f[z][1] = f_new[z][1];
            }
        }
        for (size_t c = 0; c < cluster_length; c++) {
            big_X[i][cluster[c][0]][cluster[c][1]]*=-1.0;
        }
        num_clusters++;
        clusterss += cluster_length;
    }
}

void RW_int_cluster(){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    for (double i = 1.9; i < 3.0; i += .01){
        double h_old = 0.0;
        bool gate = true;
        while (gate) {
            reset_string_cluster(i);
            double h = get_H_hat()/pow(N,2);
            if (h > h_old || h == h_old) gate = false;
            h_old = h;
            //cout << h << "\n";
        }
        for (size_t y = 0; y < 10; y++) {
            reset_string_cluster(i);
        }
        double H_hat  = get_H_hat();
        cout << i << "\t" <<  H_hat/pow(N,2) << "\t" << get_C_hat(i,H_hat) << "\t" << clusterss/num_clusters << "\t" << get_m() << "\t" << get_B4() << "\t" << get_chi(i) << "\n";
        fs << i << "\t" <<  H_hat/pow(N,2) << "\t" << get_C_hat(i,H_hat) << "\t" << clusterss/num_clusters << "\t" << get_m() << "\t" << get_B4() << "\t" << get_chi(i) << "\n";
        set_string();
        num_clusters = 0.0; clusterss = 0.0;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    set_string();
    RW_int_cluster();
    return 0;
}

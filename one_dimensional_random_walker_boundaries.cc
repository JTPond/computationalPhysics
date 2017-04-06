//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

double PI = 2. * asin(1.0);

int getRandomNumber(int min, int max){ 
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);  
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

int X(int t){
    if(t == 1) return 1;
    else return -1;
}

double pi(double x, double t){
    return sqrt(2.0/(PI*t))*exp(-(x*x)/(2.0*t));
}

double prob_after_t(double t, int L){
    int walkers[10000];
    for (int i = 0; i < 10000; i += 1)walkers[i] = 1;
    for (int j = 0; j < 10000; j += 1){
        double sum = 0.0;
        for (int k = 0; k < t; k += 1){
            double x = X(getRandomNumber(0,1));
            if(sum < L && sum > -L) sum += x;
            else{
                walkers[j] = 0;
                break;
            }
        }
    }
    double prSum = 0.0;
    for (int i = 0; i < 10000; i += 1)prSum+=walkers[i];
    return prSum/10000.0;
}

double simulate_walker_avgt(int L){
    double t =0.0;
    for (int c = 0; c < 10000; c+= 1){
        double sum = 0.0;
        for (int i = 0; i < 200000; i += 1){
            double x = X(getRandomNumber(0,1));
            if(sum < L && sum > -L) sum += x;
            else{
                t+=i;
                break;
            }
        }
    }
    return t/10000.0;
}    

void find_p_of_t(void){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/f.dat");
    for (double t = 0.0; t <= 12000.0; t += 10.0){
        fs << t << '\t';
        for (int l = 10; l <= 60; l += 10){
            fs << prob_after_t((double) t,l) << '\t';
        }
        fs << '\n';
    }
    fs.close();
}

//trapazoidal method
void RW_int(void){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/e.dat");
    double n_max = 100;
    fs.precision(8);
    //loop over possible values of N
    for(int n=0; n<=n_max; n+=10){ 
        fs << n << '\t' << simulate_walker_avgt(n) << endl;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    //RW_int();
    find_p_of_t();
    return 0;
}



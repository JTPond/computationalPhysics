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

double X(int t){
    if(t == 1) return 1.0;
    else return -1.0;
}

double pi(double x, double t){
    return sqrt(2.0/(PI*t))*exp(-(x*x)/(2.0*t));
}

//trapazoidal method
void RW_int(int T){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/c.dat");
    double n_max = 400000;
    fs.precision(8);
    //loop over possible values of N
    for(int n=1; n<=n_max; n++){ 
        double sum = 0.0;
        for (int i = 0; i < T; i += 1){
            double x = X(getRandomNumber(0,1));
            sum += x;
        }
        fs << n << '\t' << sum << endl;
    }
    fs.close();
}

double sigmax(double T){
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < T; i += 1){
        double x = X(getRandomNumber(0,1));
        sum1+=x;//*pi(x,T);
        sum2+=x*x;//*pi(x,T);
    }
    return sqrt((sum2) - (sum1*sum1));
}

void f_sig(void){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/d.dat");
    fs.precision(8);
    for (double i = 0.0; i < 5500.0; i += 50.0){
        fs << i << '\t' << sigmax(i) << endl;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    //RW_int(1000);
    f_sig();
    return 0;
}



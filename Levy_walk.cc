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

double cauchy(double x){
    return (.5/PI)/(x*x + .25);
}    

double get_random_length(){
    while(true){
        double X0 = (((double) getRandomNumber(0,25000000))/10000.);
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);
        double Y0 = cauchy(X0);
        if (X1 < Y0)return X0;
    }
}

void RW_int(){
    ofstream fs;
    fs.open("./a.dat");
    double n_max = 250000;
    fs.precision(8);
    //loop over possible values of N
    for(int n=1; n<=n_max; n++){ 
        double sum = 0.0;
        for (int i = 0; i < 5; i += 1){
            double x = get_random_length();
            sum += x;
        }
        fs << n << '\t' << sum << endl;
    }
    fs.close();
}


int main(int argc, char const *argv[]) {
    srand(time(0));
    RW_int();
    return 0;
}



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

//returns value of the function being integrated
double H(double x1, double x2, double x3, double x4){
    // No need to add 0.0, and the total integrated sum starts with 1.0 term
    if(sqrt(x1*x1 + x2*x2 + x3*x3+ x4*x4) <= 1.0) return 1.0;
    else return 0.;
}

//trapazoidal method
void MC_int(void){
    double exact_sum = (PI*PI)/2.0;

    ofstream fs;
    fs.open("/home/jpond/Desktop/workspace/ComputationalPhysics/dataFiles/a.dat");

    double n_max = 100000000;
    fs.precision(8);
    //loop over possible values of N
    double sum = 0.0;
    for(double n=1.0; n<=n_max; n++){ 
        //start with the (1/2)*f_n term        
        //Loop for actual sum
        sum += H((((double) getRandomNumber(0,20000000))/10000000.)-1.0,(((double) getRandomNumber(0,20000000))/10000000.)-1.0,(((double) getRandomNumber(0,20000000))/10000000.)-1.0,(((double) getRandomNumber(0,2000))/1000.)-1.0);
        //factor in h
        if(fmod(n,1000000) == 100000) {
            double sum1 = sum*16.0/n; 
            // write to file with abs(error)
            fs << n << '\t' << sum1 << '\t' << sqrt((sum1 - exact_sum)*(sum1 - exact_sum))/exact_sum << endl;
        }
    }

    fs.close();
}
void MC_int_rI(void){
    double exact_sum = (PI*PI)/2.0;

    ofstream fs;
    fs.open("/home/jpond/Desktop/workspace/ComputationalPhysics/dataFiles/b.dat");

    fs.precision(8);
    //loop over possible values of N
    for(double n=1; n<=150; n++){ 
        double sum = 0.0;
        for (double i = -1.0; i <= 1.0; i += 2.0/n){
            for (double j = -1.0; j <= 1.0; j += 2.0/n){
                for (double k = -1.0; k <= 1.0; k += 2.0/n){
                    for (double l = -1.0; l <= 1.0; l += 2.0/n){
                        sum+=H(i,j,k,l);
                    }
                }//endk
            }//endj
        }//endi
        double sum1 = (sum*16.0)/pow(n,4.0);
        fs << pow(n,4.0) << '\t' << sum1 << '\t' << sqrt((sum1 - exact_sum)*(sum1 - exact_sum))/exact_sum << endl;        
    }//endn
    fs.close();
}


int main(int argc, char const *argv[]) {
    //srand(time(0));
    MC_int_rI();

    return 0;
}



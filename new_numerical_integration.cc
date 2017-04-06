//C++ STD  
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

//returns value of the function being integrated
double f(double x){
    // No need to add 0.0, and the total integrated sum starts with 1.0 term
    if (x == 0.0 || x == 1.0)return 0.0;
    else return sqrt(1.0-(x*x));
}
//trapazoidal method
void trap_int(void){
    double a = 0.0;
    double b = 1.0;
    double exact_sum = PI/4.0;

    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");

    int n_max = 10000;
    fs.precision(20);
    //loop over possible values of N
    for(int n=10; n<=n_max; n+=100){ 
        //start with the (1/2)*f_n term
        double sum = .5;
        double N = double(n);
        double h = b/N;
        //Loop for actual sum
        for (int i = 0; i <= N; i++){
            double I = double(i);
            // h = 1/N, so I*h = a step between 0.0 and 1.0
            sum += I*h*f(I*I*h*h);
        }
        //factor in h
        sum *=2.0*h; 
        // write to file with abs(error)
        fs << n << '\t' << sum << '\t' << sqrt((sum - exact_sum)*(sum - exact_sum))/exact_sum << endl;
    }

    fs.close();
}
//simpson's method (redundent comments omitted)
void simp_int(void){
    double a = 0.0;
    double b = 1.0;
    double exact_sum = PI/4.0;

    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");

    int n_max = 10000;
    fs.precision(16);
    for(int n=10; n<=n_max; n+=100) {
        double sum = 1.0;
        double N = double(n);
        double h = b/N;
        for (int i = 0; i <= N; i++){
            double I = double(i);
            // if odd return 4.0*f
            if (i%2 != 0) sum += 4.0*f(I*h);
            // if even return 2.0*f
            else sum += 2.0*f(I*h);
        }
        // factor in h/3
        sum *=h/3.0; 
        fs << n << '\t' << sum << '\t' << sqrt((sum - exact_sum)*(sum - exact_sum))/exact_sum << endl;
    }

    fs.close();
}
int main (int argc, char *argv[]){
    trap_int();
    simp_int();

    return 0;
}


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

double PI = 2. * asin(1.0);

// Riemann Zeta Function
void error_analysis_s1(void) {

    double exact_sum = pow(PI, 2) / 6.;

    ofstream fs;
    fs.open("./dataFiles/a.dat");

    int n_max = 200000;
    double sum = 0;

    fs.precision(16);
    // loop for sum
    for(int n=1; n<n_max; n++) {
        // add 1/(n^2)
        sum += 1. / pow((double) n, 2);
        //write to file in 10000 n intervals 
        if(n % 10000 == 1000) {
            fs << n << '\t' << sum << '\t' << (exact_sum - sum)/exact_sum << endl;
        }
    }

    fs.close();
}
// Harmonic Series (redundent comments omitted)
void error_analysis_s2(void) {

    double exact_sum = log(2.0);

    ofstream fs;
    fs.open("./dataFiles/b.dat");

    int n_max = 200000;
    double sum = 0;

    fs.precision(16);
    for(int n=1; n<n_max; n++) {
        // add ((-1)^(n-1))/n
        sum += pow(-1.0,n-1) / n;

        if(n % 10000 == 1000) {
            fs << n << '\t' << sum << '\t' << (exact_sum - sum)/exact_sum << endl;
        }
    }

    fs.close();
}
// Geometric Series (redundent comments omitted)
void error_analysis_s3(void) {

    double r  = 1.0001;
    double exact_sum = 1.0 / (r-1.0);

    ofstream fs;
    fs.open("./dataFiles/c.dat");
    fs.precision(16);
    // do each sum indevidually, 
    for(int k=0; k<20; k++){
        double sum = 0.0;
        //do largest sum first
        int n_max = 200000 - (k*10000);
        // sum backwards 
        for(int n=n_max; n>0; n--) {
            double N = double(n);
            // add 1/r^n
            sum += 1.0 / pow(r,N);
        }
        //cout << 1.0 / pow(r, (double) n) << endl;
        // write to file for every value of k
        fs << 200000 - (k*10000) << '\t' << sum << '\t' << (exact_sum - sum)/exact_sum << endl;
    }
    fs.close();
}

int main(int argc, const char * argv[]) {

    error_analysis_s1();
    error_analysis_s2();
    error_analysis_s3();

    return 0;
}

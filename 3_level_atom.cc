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

double n(double omega,double T){
    return 1.0/(exp(omega/T)-1.0);
}    

double prob(int s, int st, double T){
    double a=0.02;
    double b=0.03;
    double g=0.01;
    double E0=0.0;
    double E1=1.0;
    double E2=2.1;
    if(s == 0 && st == 1)return a*n(E1-E0,T);
    else if(s == 0 && st == 2)return b*n(E2-E0,T);
    else if(s == 1 && st == 2)return g*n(E2-E1,T);
    else if(s == 1 && st == 0)return a*(n(E1-E0,T)+ 1.0);
    else if(s == 2 && st == 0)return b*(n(E2-E0,T)+ 1.0);
    else if(s == 2 && st == 1)return g*(n(E2-E1,T)+ 1.0);
    else if(s == 0 && st == 0)return 1.0 - a*n(E1-E0,T) - b*n(E2-E0,T);
    else if(s == 1 && st == 1)return 1.0 - g*n(E2-E1,T) - a*(n(E1-E0,T)+ 1.0);
    else if(s == 2 && st == 2)return 1.0 - b*(n(E2-E0,T)+ 1.0) - g*(n(E2-E1,T)+ 1.0);
    else return 0.0;
}

int get_random_transition(int s, double T){
    while(true){
        int X0 = getRandomNumber(0,2);
        double X1 = (((double) getRandomNumber(0,10000000))/10000000.);
        double Y0 = prob(s,X0,T);
        if (X1 < Y0)return X0;
    }
}

void RW_int(double T){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    double n_max = 250;
    fs.precision(8);
    //loop over possible values of N
    for(int n=0; n<=n_max; n++){
        double pi0 = 0.0;
        double pi1 = 0.0;
        double pi2 = 0.0;
        for(double j = 0 ; j < 10000; j++){
            int s = 0;        
            for (int i = 0; i < n; i += 1){
                s = get_random_transition(s,T);
            }
            switch(s){
                case 0:
                    pi0+=1.0;
                    break;
                case 1:
                    pi1+=1.0;
                    break;
                case 2:
                    pi2+=1.0;
                    break;
            }
        }
        fs << n << '\t' << pi0/10000.0 << '\t' << pi1/10000.0 << '\t' << pi2/10000.0 << endl;
    }
    fs.close();
}

void RW_Temp(){
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/b.dat");
    fs.precision(8);
    for(double T=0.0; T<=2.0; T+=.125){
        double pi0 = 0.0;
        double pi1 = 0.0;
        double pi2 = 0.0;
        for(double j = 0 ; j < 10000; j++){
            int s = 0;        
            for (int i = 0; i < 1500; i += 1){
                s = get_random_transition(s,T);
            }
            switch(s){
                case 0:
                    pi0+=1.0;
                    break;
                case 1:
                    pi1+=1.0;
                    break;
                case 2:
                    pi2+=1.0;
                    break;
            }
        }
        fs << T << '\t' << (pi1/10000.0)/(pi0/10000.0) << '\t' << (pi2/10000.0)/(pi0/10000.0) << endl;
    }
    fs.close();
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    RW_int(1.4);
    RW_Temp();
    return 0;
}



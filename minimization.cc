//C++ STD
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

int getRandomNumber(int min, int max){
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

double x0,x1,x2;
double z0,z1,z2;
double xh,zh;

double function(double x, double z){
    double a = 1;
    double b = 10;
    return pow(x-a,2.0) + b*pow(z - pow(x,2.0),2.0);
}

void order(){
    double X0 = function(x0,z0);
    double X1 = function(x1,z1);
    double X2 = function(x2,z2);
    if(X2 > X1 && X2 > X0){
        xh=x0; zh=z0;
        x0=x2; z0=z2;
        x2=xh; z2=zh;
        if(X2 > X1){
            xh=x1; zh=z1;
            x1=x2; z1=z2;
            x2=xh; z2=zh;
        }
    }
    else if(X1 > X2 && X1 > X0){
        xh=x0; zh=z0;
        x0=x1; z0=z1;
        x1=xh; z1=zh;
        if(X2 > X1){
            xh=x1; zh=z1;
            x1=x2; z1=z2;
            x2=xh; z2=zh;
        }
    }
    else if(X0 > X2 && X0 > X1){
        if(X2 > X1){
            xh=x1; zh=z1;
            x1=x2; z1=z2;
            x2=xh; z2=zh;
        }
    }
}

void cycle(){
    double _x = (x2 + x1)/2.0;
    double _z = (z2 + z1)/2.0;
    double _vx = (_x-x0);
    double _vz = (_z-z0);
    double xa = x0 + 2*_vx; double za = z0 + 2*_vz;
    double xb = x0 + 3*_vx; double zb = z0 + 3*_vz;
    double xc = x0 + 1.5*_vx; double zc = z0 + 1.5*_vz;
    double xd = x0 + .5*_vx; double zd = z0 + .5*_vz;
    double XA = function(xa,za);
    if(XA < function(x2,z2)){
        double XB = function(xb,zb);
        if(XA < XB){
            x0=xa; z0=za;
        }
        else if (XB < XA){
            x0=xb; z0=zb;
        }
    }
    else if(XA > function(x2,z2) && XA < function(x1,z1)){
        x0=xa; z0=za;
    }
    else if(XA > function(x1,z1)){
        double XC = function(xc,zc);
        if(XC < function(x1,z1)){
            x0=xc; z0=zc;
        }
        else if (XC > function(x1,z1)){
            double XD = function(xd,zd);
            if(XD < function(x1,z1)){
                x0=xd; z0=zd;
            }
            else if (XD > function(x1,z1)){
                x0 = (x2+x0)/2.0; z0 = (z2+z0)/2.0;
                x1 = (x2+x1)/2.0; z1 = (z2+z1)/2.0;
            }
        }
    }
}

int main(int argc, char const *argv[]) {
    srand(time(0));
    x0 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    x1 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    x2 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    z0 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    z1 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    z2 = (((double) getRandomNumber(0,10000000))/1000000. - 5.0);
    order();
    double __x = (x2 + x1 + x0)/3.0;
    double __z = (z2 + z1 + z0)/3.0;
    double __f = function(__x,__z);
    double __nf;
    cout << __x << ','<< __z << ','<< __f << ',' << '\n';
    int i = 0;
    int r = 0;
    ofstream fs;
    fs.open("/home/jpond/workspace/ComputationalPhysics/dataFiles/a.dat");
    while(r<50){
        cycle();
        order();
        __x = (x2 + x1 + x0)/3.0;
        __z = (z2 + z1 + z0)/3.0;
        __nf = function(__x,__z);
        fs << i << ' ' << __x << ' '<< __z << ' '<< __f << ' ' << '\n';
        if(__nf == __f){
            r++;
        }
        __f = __nf;
        i++;
    }
    fs.close();

    return 0;
}

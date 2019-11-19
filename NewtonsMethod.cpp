#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

const double threshhold = 0.001;
const int maxits = 1000;

double funct(double x){
    return pow(x,3) - 4*pow(x,2) + x + 3;
}
double functdx(double x){
    return 3*pow(x,2) - 8*x + 1;
}

double newtons(double (*f)(double), double (*dfdx)(double), double x0){
    double xn = x0;
    double fxn,dfxn;
    printf("  n  |    xn    |   f(xn)  | dfdx(xn) \n--------------------------------------\n");

    for(int n = 0; n < maxits; n++){
        fxn = (*f)(xn);
        dfxn = (*dfdx)(xn);
        printf(" %3d | %8.4f | %8.4f | %8.4f \n", n, xn, fxn, dfxn);
        if(abs(fxn) < threshhold){
            printf("f(%f) = %f\n", xn, fxn);
            return xn;
        }
        xn = xn - (fxn / dfxn);
    }
    cout << "did not converge\n";
    return 0;
}

int main(int argc, char* argv[]){
    double x0 = 5;
    newtons(funct, functdx, x0);

    return 0;
}

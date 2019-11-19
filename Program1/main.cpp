#include <math.h>
#include <iostream>

using namespace std;

/*
Symmetric Derivative
= lim(h -> 0) of [f(x+h)-f(x-h)]/2h
= [f(x + delta_t) - f(x - delta_t)]/(d*delta_t)

*/

double Derivative(double (*f)(double), const double delta_t, const double x){
    return ((*f)(x + delta_t) - (*f)(x - delta_t))/(2*delta_t);
}



double absvalue(double x){
    return ((x < 0)?(-1*x):x);
}
double derivabsvalue(double x){
    return ((x!=0)?((x < 0)?-1:1):0);
}

double parabola(double x){
    return (x*x);
}
double derivparabola(double x){
    return 2*x;
}

double polynomial(double x){
    return 2*x*x*x - 3*x*x + 0.4*x - 0.7;
}
double derivpolynomial(double x){
    return 6*x*x - 6*x + 0.4;
}

int main(int argc, char* argv[]){
    double (* functs [])(double) = {absvalue,parabola,polynomial};
    double (* derivs [])(double) = {derivabsvalue,derivparabola,derivpolynomial};
    const string fnames[] = {"|x|","x^2","2x^3 - 3x^2 + 0.4x - 0.7"};
    const unsigned int num_fs = 3;
    const double x_vals[] = {-5,-1,0,1,5};
    const unsigned int num_xs = 5;
    const double t_vals[] = {1,0.1,0.01,0.001,0.0001};
    const unsigned int num_ts = 5;

    double actual,approx;

    for(int p = 0; p < num_xs; p++){
        cout << "x = " << x_vals[p] << "\n\n";
        for(int n = 0; n < num_fs; n++){
            actual = (*(derivs[n]))(x_vals[p]);
            cout << "\tf(x) = " << fnames[n] << "  actual:" << actual <<"\n";
            printf("\t\t delta_t | approx   | error\n");
            for(int i = 0; i < num_ts; i++){
                approx = Derivative(functs[n],t_vals[i],x_vals[p]);
                printf("\t\t%6f | %4.4f | %4.4f\n",t_vals[i],approx,abs(actual-approx));
            }
        }
    }

    return 0;
}

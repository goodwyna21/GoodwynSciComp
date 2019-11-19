#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;
const string DEF_FNAME = "Program2Out.txt";

//these methods take a function pointer f
double leftEndpoint(double a, double b, int n, double (*f)(double)){
    double sum = 0;
    for(double i = a; i < b; i+=(b-a)/n){
        sum+=(*f)(i);
    }
    return ((b-a)/n)*sum;
}

double rightEndpoint(double a, double b, int n, double (*f)(double)){
    double sum = 0;
    for(double i = a+((b-a)/n); i <= b; i+=((b-a)/n)){
        sum+=(*f)(i);
    }
    return  ((b-a)/n)*sum;
}

double trapezoidal(double a, double b, int n, double (*f)(double)){
    double sum = (*f)(a) + (*f)(b);
    for(double i = a+((b-a)/n); i<b; i+=((b-a)/n)){
        sum+=2*(*f)(i);
    }
    return ((b-a)/(2*n))*sum;
}

//actual functions
double funct1(double x){
    return x;
}

double funct2(double x){
    return x*x;
}

double funct3(double x){
    return x*x*x;
}

double funct4(double x){
    return (3*x*x)-(4*x)+3;
}

//hand-computed integrals
double integral1(double a, double b){
    return ((b*b)-(a*a))/2;
}
double integral2(double a, double b){
    return ((b*b*b)-(a*a*a))/3;
}
double integral3(double a, double b){
    return ((b*b*b*b)-(a*a*a*a))/4;
}

int main(int argc, char* argv[]){
    double a = 0;
    double b = 4;
    ofstream outf((argc>=2)?argv[1]:DEF_FNAME);
    printf("Integral from a=%.3f to b=%.3f\n\n",a,b);

    /*
    For the purposes of displaying, I used lots of arrays to avoid rewriting code
    */
    double (* functions [])(double) = {funct1,funct2,funct3};
    double (* integrals [])(double,double) = {integral1,integral2,integral3};
    double (* methodFuncts [])(double,double,int,double (*f)(double)) = {leftEndpoint,rightEndpoint,trapezoidal};
    int steps[] = {1,5,10,15,50};
    int numSteps = 5;
    string methodNames[] = {"Left Endpoint","Right Endpoint","Trapezoidal"};
    string labels[] = {"f(x)=x  ","f(x)=x^2","f(x)=x^3"};
    double tmp1, tmp2;
    outf << ",0";
    for(int k = 0; k < 3; k++){ //k loops through the three methods
        printf("%s\n",methodNames[k].c_str());
        printf("Steps:     actual");
        for(int h = 0; h < numSteps; h++){ //h loops through to create the header
            printf("          %.3s      ",to_string(steps[h]).c_str());
            if(k == 0){
                outf << "," << steps[h];
            }
        }
        outf << "\n" << methodNames[k] << "\n";
        printf("\n");
        for(int i = 0; i < 3; i++){ //i loops through the three functions
            printf("%s ",labels[i].c_str());
            outf << labels[i] << ",";
            tmp1 = (*integrals[i])(a,b);
            printf("| %.6s ",to_string(tmp1).c_str()); //output actual integral
            outf << tmp1;
            for(int h = 0; h < numSteps; h++){ //h loops through the values for steps
                tmp2 = (*methodFuncts[k])(a,b,steps[h],functions[i]);
                printf("| %.6s (%.6s) ",to_string(tmp2).c_str(),to_string(tmp2-tmp1).c_str());
                outf << "," << tmp2;
            }
            printf("\n");
            outf << "\n";
        }
        printf("\n");
        outf << "\n";
    }
    outf.close();

    return 0;
}

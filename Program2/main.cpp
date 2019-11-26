#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;

//CONSTANTS
const string DEFAULT_FNAME = "Program2Out.txt";
const double A = 0;
const double B = 4;
const unsigned int NUMSTEPS = 5;
const unsigned int STEPS[NUMSTEPS] = {1,5,10,15,50};
const string METHODNAMES[] = {"Left Endpoint","Right Endpoint","Trapezoidal"};
const string LABELS[] = {"x       ","x^2     ","8-0.5x^2"}; //all 8 chars long for outputting

//testing functions
//x
double funct1(double x){
    return x;
}

//x^2
double funct2(double x){
    return x*x;
}

//-0.5(x^2) + 8
double funct3(double x){
    return -0.5*(x*x) + 8;
}

//hand-computed integrals
double actualintegral1(double a, double b){
    return ((b*b)-(a*a))/2;
}
double actualintegral2(double a, double b){
    return ((b*b*b)-(a*a*a))/3;
}
double actualintegral3(double a, double b){
    return ((1.0/12)*(b*b*b - a*a*a)) + 4*(b - a);
}
//END CONSTANTS

//these methods take a function pointer f and integrate it from a to b with the corresponding method
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

int main(int argc, char* argv[]){
    ofstream outf((argc>=2)?argv[1]:DEFAULT_FNAME); //optional command line argument specifies output file
    printf("Integral from a=%.3f to b=%.3f\n\n",A,B);

    /*
    For the purposes of displaying, I used lots of arrays to avoid rewriting code
    The output file is formatted as a .csv so that it can be imported into a spreadsheet
    There is also command line output with the same data
    */
    double (* functions [])(double) = {funct1,funct2,funct3};
    double (* integrals [])(double,double) = {actualintegral1,actualintegral2,actualintegral3};
    double (* methodFuncts [])(double,double,int,double (*f)(double)) = {leftEndpoint,rightEndpoint,trapezoidal};
    double tmp1, tmp2;

    outf << A << "," << B << "\n,0"; //This is just here to format the output file nicely

    //This section is pretty hard to read but it is used to loop through the functions and different steps and stuff to output it all nicely without rewriting lots of code
    for(int k = 0; k < 3; k++){ //k loops through the three methods
        printf("%s\n",METHODNAMES[k].c_str());
        printf("Steps:     actual");
        for(int h = 0; h < NUMSTEPS; h++){ //h loops through to create the header
            printf("          %.3s      ",to_string(STEPS[h]).c_str());
            if(k == 0){
                outf << "," << STEPS[h];
            }
        }
        outf << "\n" << METHODNAMES[k] << "\n";
        printf("\n");

        for(int i = 0; i < 3; i++){ //i loops through the three functions
            printf("%s ",LABELS[i].c_str());
            outf << LABELS[i] << ",";
            tmp1 = (*integrals[i])(A,B);
            printf("| %.6s ",to_string(tmp1).c_str()); //output actual integral
            outf << tmp1;
            for(int h = 0; h < NUMSTEPS; h++){ //h loops through the values for steps
                tmp2 = (*methodFuncts[k])(A,B,STEPS[h],functions[i]);
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

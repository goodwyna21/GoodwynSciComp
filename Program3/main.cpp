#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>

using namespace std;

const double G = 9.81;
const double PI = 3.14159265;
const double H = 0.001;

/*
Runge Kutta to integrate the equations for motion of a pendulum
dx/dt = v
dv/dt = -(G/l)sin(x)
x: angular displacement from vertical
v: angular velocity
l: length of pendulum

x(t) = integral[0,t](v(T)dT)
v(t) = integral[0,t](-(G/l)sin(x(T))dT)
x(0) = InitialX
v(0) = 0

Runge Kutta:
y(x+h)=y(x)+ 1/6(k1+2k2+2k3+k4)
k1: h*f(x,y)
k2: h*f(x+h/2,y+k1/2)
k3: h*f(x+h/2,y+k2/2)
k4: h*f(x+h,y+k3)
*/

/*struct pendulum;
double fx(pendulum* p, double vn, double xn);
double fv(pendulum* p, double xn, double vn);
double rungekutta(double (*f)(pendulum*,double, double), pendulum* p, double x, double y, double h);
*/

struct pendulum{
    double L; //length of pendulum
    double xn;
    double vn;

    pendulum(double X0, double V0, double l)
    : L(l),xn(X0),vn(V0) {}

    double updatex(double x, double v, double h){
        double k1 = h * vn;
        double k2 = h * (vn + (h/2));
        double k3 = h * (vn + (h/2));
        double k4 = h * (vn + h);
        xn = x + (1/6)*(k1 + (2*k2) + (2*k3) + k4);
    }

    double dvdt(double x){
        return -1*(G/L)*sin(x);
    }

    double updatev(double x, double v, double h){
        double k1 = h * dvdt(x);
        double k2 = h * dvdt(x + (h/2));
        double k3 = h * dvdt(x + (h/2));
        double k4 = h * dvdt(x + h);
        vn = v + (1/6)*(k1 + (2*k2) + (2*k3) + k4);
    }

    void mainloop(int steps, double h, ofstream * outf, bool verbose=false){
        printf("   step |       x |      v \n %6d | %.6f | %.6f \n",0,xn,vn);
        (*outf) << "n,x,v\n";
        double oldx,k1,k2,k3,k4;
        for(int i = 0; i < steps; i++){
            oldx = xn;

            /*k1 = h * vn;
            k2 = h * (vn + (h/2));
            k3 = h * (vn + (h/2));
            k4 = h * (vn + h);
            xn += (1.0/6)*(k1 + (2*k2) + (2*k3) + k4);*/
            xn += h*vn;
            vn += h*dvdt(oldx);
            /*
            k1 = h * dvdt(oldx);
            k2 = h * dvdt(oldx + (h/2));
            k3 = h * dvdt(oldx + (h/2));
            k4 = h * dvdt(oldx + h);
            vn += (1.0/6)*(k1 + (2*k2) + (2*k3) + k4);*/

            if(verbose){
                (*outf) << i+1 << "," << xn << "," << vn << "\n";
                printf(" %6d | %.6f | %.6f \n",i+1,xn,vn);
            }
        }
    }
};

/*
double rungekutta(double (*f)(pendulum*,double, double), pendulum* p, double a, double b, double h){
    //solves for b(x(n+1))
    double k1 = h * (*f)(p,a,b);
    double k2 = h * (*f)(p,a+(h/2),b+(k1/2));
    double k3 = h * (*f)(p,a+(h/2),b+(k2/2));
    double k4 = h * (*f)(p,a+h,b+k3);
    return (b + (1/6)*(k1 + 2*k2 + 2*k3 + k4));
}

double fx(pendulum* p, double vn, double xn){
    return vn;
}

double fv(pendulum* p, double xn, double vn){
    return -(G/(p->L))*sin(xn);
}
*/

int main(int argc, char* argv[]){
    double x0 = 0.1;
    double v0 = 0;
    double l = 1;
    int steps = 10000;
    pendulum p(x0,v0,l);
    if(argc < 2){
        cout << "Error: Need filename\n";
        return 1;
    }
    ofstream outf(argv[1], ios::out);
    printf("x0: %f\nv0: %f\n L: %f\n n: %d\n h: %f\n", x0,v0,l,steps,H);
    p.mainloop(steps,H,&outf,true);

    return 0;
}

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;


void syeuler(double* p, double* q, double dt, double H, const double dim, double pq);


int main (){
    
const int dim = 2;
double q[2], p[2];
double H=0;
const double dt = 0.0005;
const double PI = 3.14; 
const double e = 0.6;
const double tEnd = 20*PI;
const double t0 = 0;
const int N = (tEnd-t0)/dt;
double t = t0;
double pq = 0;
q[0] = 1-e;
q[1] = 0;
p[0] = 0;
p[1] = sqrt((1+e)/(1-e));

ofstream out("data.txt");
H = (p[0] * p[0] + p[1]*p[1])/2 - (1/sqrt(pow(q[0],2)+pow(q[1],2)));




out << t0 << "\t" <<  q[0] << "\t" << q[1] <<"\t"<< p[0]<< "\t" << p[1] << "\t" << H << endl;


    for (int i = 0; i<N ; i++){
    
    syeuler(q, p , dt, H , dim, pq);
    
    t += dt;
    
    out << t << "\t" <<  q[0] << "\t" << q[1] <<"\t"<< p[0]<< "\t" << p[1] << "\t" << H << endl;
    }
    out.close();
    
    return 0;
    
    
    
}
    
    void syeuler(double* p, double* q, double dt, double H, const double dim,double pq){
   
    pq=q[0]*q[0] + q[1]*q[1]; 
    for(int j=0; j < dim; j++){
        
        p[j] = p[j] - dt *q[j] * pow(pq,-3.0/2.0);
        q[j] = q[j] + dt * p[j];
        
    }
    H = (p[0] * p[0] + p[1]*p[1])/2.0 - (1/sqrt(pow(q[0],2.0)+pow(q[1],2.0)));
}
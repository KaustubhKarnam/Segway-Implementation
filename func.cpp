/* Introduction to Scientific Programming, Programming project, 2019/2020 */
/* (c) Matthias Braendel, 2019 */
/* func.cpp */

#include <iostream>
#include <fstream>
#include <vector>
#include "func.h"

void f_introduction(){
    std::cout << "\nHello, this is the programming project for the students of the module";
    std::cout << "\nIntroduction to Scientific Programming, WS 2019/2020, TUBAF";
    std::cout << "\n\nGoal is to complete the code given to create a little movie of ";
    std::cout << "an inverse pendulum balancing itself!\n";
}

// https://www.geeksforgeeks.org/passing-vector-function-cpp/
void f_vecinit(std::vector<double> &v){
    for(unsigned int m = 0; m < v.size(); m++){
        v[m] = 0.0;
    }
}

void f_filltspan(std::vector<double> &v, double tsta, double tend, int framerate){
    // #################
    // Your code here...
    // TASK 2.1: fill the tspan array with time-steps: tstart, tstart+h,...tend

    for (unsigned int i = 0; i < v.size(); i++){
        v[i] = tsta+(i*(1.0/framerate));
    }
    // #################
}

void f_printtimeoptions(const std::vector<double> &v, double tsta, double tend,int framerate){
    // #################
    // Your code here...
    // TASK 2.2: print time options as shown below, but with your actual values from your actual variables
    // #################
    double n;
    n = ((tend-tsta)*framerate)+1;
    std::cout << "\n  Expected values:\n\t start time: "<< tsta <<"\n\t end time: "<< tend <<"\n\t framerate:" << framerate <<"\n\t #steps:" << n <<"\n";
}

void f_printmodelparameters(double g, double m1, double m2, double l){
    // #################
    // Your code here...
    // TASK 2.3: see f_printtimeoptions
    // #################
    std::cout << "\n  Expected values:\n\t gravitational constant:" << g << "\n\t mass of wagon:" << m1 << "\n\t mass of pendulum:" << m2 << "\n\t length of pendulum: "<< l <<"\n";
}

void f_printcoefficients(const std::vector<double> &C){
    // #################
    // Your code here...
    // TASK 2.4: see f_printmodelparameters
    // #################
    std::cout << "\n  Expected values:\n\t" << C[0]<< "\n\t" << C[1] << "\n\t" << C[2] << "\n\t" << C[3] << "\n\t" << C[4] << "\n\t" << C[5] <<"\n";
}

void f_filldisturbance(std::vector<double> &T, double valT, std::vector<double> tdist, int framerate){
    int index = int(tdist[0]*framerate);
    for(unsigned int m = 0; m < tdist.size(); m++){
        T[index+m] = valT;
    }
}

void f_printdisturbance(double valT, double tdiststa, double tdistend, int framerate){
    // #################
    // Your code here...
    // TASK 2.5: see f_printcoefficients

    std::cout << "\n  Expected values:\n\t start time:" << tdiststa << "\n\t end time: " << tdistend << "\n\t stepsize: " << 1.0/(framerate) << "\n\t magnitude: " << valT << "\n";
    // #################
}

void f_matinit(std::vector< std::vector<double> > &A){
    unsigned int nofrows = A.size();
    for(unsigned int mr = 0; mr < nofrows; mr++){
        unsigned int nofcols = A[mr].size();
        for(unsigned int mc = 0; mc < nofcols; mc++){
            A[mr][mc] = 0.0;
        }
    }
    /* or
    for(unsigned int mr = 0; mr < nofrows; mr++){
        f_vecinit(A[mr]);
    } */
}

void f_printlinearsystem(const std::vector< std::vector<double> > &A,const std::vector<double> &b1, const std::vector<double> &b2, const std::vector<double> &K){
    std::cout << "\nLinear system A*x+b1*F+b2*T:";
    // #################
    // Your code here...
    // TASK 2.6: obtain the following format with your variables
    std::cout << "\n  Expected values:\n\t|0.000	1.000	0.000	0.000||x_1|  +  |0.000|*F  +  |0.000|*T\n\t|0.000	0.000	2.264	0.000||x_2|  +  |0.855|*F  +  |3.142|*T\n\t|0.000	0.000	0.000	1.000||x_3|  +  |0.000|*F  +  |0.000|*T\n\t|0.000	0.000	44.389	0.000||x_4|  +  |3.142|*F  +  |61.613|*T\n";
    std::cout << "\n Generated Values: \n\t";
    std::cout << "|";
    for (unsigned int i = 0; i < A[0].size(); i++) {
        std::cout << A[0][i] << "\t";
    }
    std::cout << "\t";
    std::cout << "| |x_1| + ";
    std::cout << "|" << b1[0] << "|*F \t+ \t" << "|" << b2[0] << "|*T \n\t";

    std::cout << "|";
    for (unsigned int i = 0; i < A[1].size(); i++) {
        std::cout << A[1][i] << "\t";
    }
    std::cout << "| |x_2| + ";
    std::cout << "|" << b1[1] << "|*F + " << "|" << b2[1] << "|*T \n\t";

    std::cout << "|";
    for (unsigned int i = 0; i < A[2].size(); i++) {
        std::cout << A[2][i] << "\t";
    }
    std::cout << "\t";
    std::cout << "| |x_3| + ";
    std::cout << "|" << b1[2] << "|*F \t+ \t" << "|" << b2[2] << "|*T \n\t";

    std::cout << "|";
    for (unsigned int i = 0; i < A[3].size(); i++) {
        std::cout << A[3][i] << "\t";
    }
    std::cout << "| |x_4| + ";
    std::cout << "|" << b1[3] << "|*F + " << "|" << b2[3] << "|*T \n\t";
    // #################
    std::cout << "\nControl Matrix K:\n\t|";
    // #################
    // Your code here...
    // TASK 2.6: display the control matrix K, that is actually a vector
    // #################
    std::cout << "|\n  Expected values:\n\t|-0.779 -1.622 25.477 3.624|\n";
    std::cout << "|\n  Generated values:\n\t|" << K[0] <<"\t"<< K[1] <<"\t"<<K[2] <<"\t"<<K[3] <<"|\n";

}

void f_systemhasbeensetup(){
    std::cout << "\nThe system has been set up successfully!\n";
}

unsigned int f_checksize(const unsigned int v, const unsigned int w){
    /* check for equal size of the vectors, if not ... */
    if(v != w){
        std::cout << "Different vector sizes!\n";
        /* check, which vector is smaller, take its size ... */
        if(v < w){
            return v;
        }else{
            return w;
        }
    }else{
        return v;
    }
}

void f_veccpy(std::vector<double> &v, const std::vector<double> &w){
    unsigned int msize = f_checksize(v.size(),w.size());
    // #################
    // Your code here...
    // TASK 5.1: copy v's values on w
    for(unsigned int i = 0; i < msize; i++){
        v[i] = w[i];
    }
    // #################
}

double f_scalarp(const std::vector<double> &v, const std::vector<double> &w){
    unsigned int msize = f_checksize(v.size(),w.size());
    double val = 0.0;
    // #################
    // Your code here...
    // TASK 5.2: do the scalar product computation
    for (unsigned i = 0; i < msize; i++){
        val += v[i]*w[i];
    }
    // #################
    return val;
}

void f_matvecmult(std::vector<double> &w, const std::vector< std::vector<double> > &A, const std::vector<double> &v){
    /* row size: */
    unsigned int mrsize = f_checksize(A.size(),w.size());
    /* col size: */
    unsigned int mcsize = f_checksize(A[0].size(),v.size());
    /* w = A*v */
    // #################
    // Your code here...
    // TASK 5.3: do the matrix vector multiplication
    for (unsigned int i = 0; i < mrsize; i++){
        for (unsigned int j = 0; j < mcsize; j++){
            w[i] += v[j]*(A[i][j]);
        }
    }
    // #################
}

void f_vecupd(std::vector<double> &w, const double a, const std::vector<double> &v){
    // #################
    // Your code here...
    // TASK 5.4: write a routine, that assigns the vector w = w + a*v
    unsigned int msize = f_checksize(v.size(), w.size());
    f_vecinit(w);
    for (unsigned int i = 0; i < msize; i++){
        w[i] = w[i]+a*v[i];
    }
    // #################
}

void f_succeed(){
    std::cout << "\nThe classical Runga-Kutta method finished!\n";
}

void f_fileoutput(const std::vector<double> &tspan, const std::vector< std::vector<double> > &y, const std::vector<double> &F, const std::vector<double> &T){
    std::ofstream F_out;
    F_out.open("segway.txt", std::ios::out);
    for(unsigned int m = 0; m < tspan.size(); m++){
        F_out << tspan[m] << " " << y[0][m+1] << " " << y[2][m+1] << " " << F[m] << " " << T[m] << "\n";
    }
    F_out.close();
    std::cout << "\nWrote data to file segway.txt!\n";
}

void f_vecshow(const std::vector<double> &v){
    for(unsigned int m = 0; m < v.size(); m++){
        std::cout << v[m] << "\n";
    }
}

void f_matshow(const std::vector< std::vector<double> > &A){
    for(unsigned int i = 0; i < A.size(); i++){
        for(unsigned int j = 0; j < A[i].size(); j++){
            std::cout << A[i][j] << " ";
        }
    std::cout << "\n";
}
}

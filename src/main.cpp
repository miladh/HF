#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;



double kineticIntegral(const double a, const double b, rowvec R_1, rowvec R_2);
double nuclearAttarctionIntegral(const double a, const double b, const int Z, rowvec R_1, rowvec R_2);
double overlapIntegral(const double a, const double b, rowvec R_1, rowvec R_2);
double errorFunction(double arg);

int main()
{
    /*------------------------------------------------------------*/
    //Constants:
//    double alpha_1 = 0.297104;
//    double alpha_2 = 1.236745;
//    double alpha_3 = 5.749982;
//    double alpha_4 = 38.216677;

    double alpha_1 = 13.00773;
    double alpha_2 = 1.962079;
    double alpha_3 = 0.444529;
    double alpha_4 = 0.1219492;

    rowvec alphas = zeros(1,4);
    alphas(0) = alpha_1;
    alphas(1) = alpha_2;
    alphas(2) = alpha_3;
    alphas(3) = alpha_4;

    int Z = 1.0;
    int numGTO = 8;

    mat h = zeros(numGTO,numGTO);
    mat S = zeros(numGTO,numGTO);
    mat F = zeros(numGTO,numGTO);

    rowvec3 R_1 = zeros(1,3);
    rowvec3 R_2 = zeros(1,3);
    mat Rij = zeros(4,3);
    R_2(0) = 1.0;

    Rij.row(1) = R_1 - R_2;
    Rij.row(2) = R_2 - R_1;


    /*------------------------------------------------------------*/


    //Set up the h matrix:
    for(uint i=0; i <alphas.size(); i++){
        for(uint j=0; j <alphas.size(); j++){

            h(i,j) = kineticIntegral(alphas[i],alphas[j],R_1,R_1)
                    + nuclearAttarctionIntegral(alphas[i],alphas[j],Z,R_1,R_1);

            h(i+4,j) = kineticIntegral(alphas[i],alphas[j],R_2,R_1)
                    + nuclearAttarctionIntegral(alphas[i],alphas[j],Z,R_2,R_1);

            h(i,j+4) = kineticIntegral(alphas[i],alphas[j],R_1,R_2)
                    + nuclearAttarctionIntegral(alphas[i],alphas[j],Z,R_1,R_2);

            h(i+4,j+4) = kineticIntegral(alphas[i],alphas[j],R_2,R_2)
                    + nuclearAttarctionIntegral(alphas[i],alphas[j],Z,R_2,R_2);


            S(i,j) = overlapIntegral(alphas[i],alphas[j],R_1,R_1);
            S(i+4,j+4) = S(i,j);

            S(i+4,j) = overlapIntegral(alphas[i],alphas[j],R_2,R_1);
            S(i,j+4) = S(i+4,j);

        }
    }

    /*------------------------------------------------------------*/
    //Set up the F matrix
    F = h;

    /*------------------------------------------------------------*/
    //Diagonalize S:

    vec s;
    mat U;
    eig_sym(s, U, S);

    mat V = U*diagmat(1.0/sqrt(s));

    F = V.t() * F * V;


    vec E;
    mat Evec;
    eig_sym(E, Evec, F);

    /*------------------------------------------------------------*/



    cout << E << endl;
    cout <<"Energy: " << E.min() << endl;

    return 0;
}


/*------------------------------------------------------------*/
double errorFunction(double arg){

    if (arg <= 1e-6){
        return 1.0;
    }

    else{
        arg = sqrt(arg);
        double f = 1.0/arg * erf(arg) *sqrt(acos(-1))/2.0;
        return f;
    }

}


/*------------------------------------------------------------*/
double nuclearAttarctionIntegral(const double a, const double b,const int Z, rowvec R_1,rowvec R_2){

    double factor = 1.0/(a+b);
    double R12 =dot(R_1-R_2,R_1-R_2);
    double expTerm = exp(-a*b*factor*R12);

    rowvec Rp = (a*R_1 + b*R_2)*factor;

    mat Rc = zeros(2,3);
    Rc(1,0) = 1.0;

    double F0_1 = errorFunction(1.0/factor*dot(Rp-Rc.row(0),Rp-Rc.row(0)));
    double F0_2 = errorFunction(1.0/factor*dot(Rp-Rc.row(1),Rp-Rc.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0_1+F0_2);

    return nucAtt;

}
/*------------------------------------------------------------*/
double kineticIntegral(const double a, const double b, rowvec R_1, rowvec R_2){

    double factor = a*b/(a+b);
    double R12 =dot(R_1-R_2,R_1-R_2);
    double expTerm = exp(-factor*R12);
    double kin = 0.5*factor*(6-4*factor*R12)
            *pow(acos(-1)/(a+b),3.0/2.0)*expTerm;

    return kin;
}


/*------------------------------------------------------------*/
double overlapIntegral(const double a, const double b, rowvec R_1, rowvec R_2){

    double factor = 1.0/(a+b);
    double R12 =dot(R_1-R_2,R_1-R_2);
    double expTerm = exp(-a*b*factor*R12);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;

}



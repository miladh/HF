#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;


double nuclearRepulsion();
double kineticIntegral(const double p, const double q, rowvec Rp, rowvec Rq);
double nuclearAttarctionIntegral(const double p, const double q, const int Z, const int Rp, const int Rq, const mat R);
double overlapIntegral(const double p, const double q, rowvec Rp, rowvec Rq);
double errorFunction(double arg);
vec normalize(vec C);
double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   rowvec Rp, rowvec Rr, rowvec Rq, rowvec Rs,
                                   vec alpha, mat S);

int main()
{

    /*-----------------------------------------------------------------------------------------------------------*/


    //System configuration:

    uint nBasisFunc = 4;
    uint nNuclei    = 2;
    uint nElectrons = 1;
    uint nSteps     = 20;

    uint Z = 1;

    vec alpha = zeros(nBasisFunc);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;

    //Initilize:
    uint nOrbitals = nBasisFunc * nNuclei;

    mat h = zeros(nOrbitals,nOrbitals);
    mat G = zeros(nOrbitals,nOrbitals);
    mat S = zeros(nOrbitals,nOrbitals);
    mat F = zeros(nOrbitals,nOrbitals);
    mat R = zeros(nNuclei,3);
    vec C = ones(nOrbitals)*0.125;

    double ****Q;
    Q = new double***[nOrbitals];
    for (uint i = 0; i < nOrbitals; ++i) {
        Q[i] = new double**[nOrbitals];

        for (uint j = 0; j < nOrbitals; ++j){
            Q[i][j] = new double*[nOrbitals];

            for (uint k = 0; k < nOrbitals; ++k){
                Q[i][j][k] = new double[nOrbitals];
            }
        }
    }

    //One nuclei at origo, the other a distance d apart, along the x-axis
    R(1,0) = 1.0;

    /*-----------------------------------------------------------------------------------------------------------*/

    //Set up the h and S matrix:
    for(uint Rp = 0; Rp < R.n_rows; Rp++){
        for(uint Rq = 0; Rq < R.n_rows; Rq++){

            for(uint p=0; p <alpha.size(); p++){
                for(uint q=0; q <alpha.size(); q++){

                    S(p+Rp*4,q+Rq*4) = overlapIntegral(alpha[p],alpha[q],R.row(Rp),R.row(Rq));
                    h(p+Rp*4,q+Rq*4) = kineticIntegral(alpha[p],alpha[q],R.row(Rp),R.row(Rq))
                            + nuclearAttarctionIntegral(alpha[p],alpha[q],Z,Rp,Rq,R);
                }
            }
        }

    }

    /*-----------------------------------------------------------------------------------------------------------*/

    //Set up the Q array:
    for(uint Rp = 0; Rp < R.n_rows; Rp++){
        for(uint Rr = 0; Rr < R.n_rows; Rr++){
            for(uint Rq = 0; Rq < R.n_rows; Rq++){
                for(uint Rs = 0; Rs < R.n_rows; Rs++){

                    for(uint p=0; p <alpha.size(); p++){
                        for(uint r=0; r <alpha.size(); r++){
                            for(uint q=0; q <alpha.size(); q++){
                                for(uint s=0; s <alpha.size(); s++){

                                    Q[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral(p,r,q,s,R.row(Rp),R.row(Rr),R.row(Rq),R.row(Rs),alpha,S);
                                }
                            }
                        }
                    }

                }
            }
        }
    }


    /*-----------------------------------------------------------------------------------------------------------*/


    for(uint step=0; step < nSteps; step++){
        G = zeros(nOrbitals,nOrbitals);

        //Set up the G matrix:
        for(uint p=0; p < nOrbitals; p++){
            for(uint q=0; q < nOrbitals; q++){
                for(uint r=0; r < nOrbitals; r++){
                    for(uint s=0; s < nOrbitals; s++){
                        G(p,q) += Q[p][r][q][s]*C(r)*C(s);
                    }
                }
            }
        }

        /*-----------------------------------------------------------------------------------------------------------*/
        //Set up the F matrix
        F = h + G;

        /*-----------------------------------------------------------------------------------------------------------*/
        //Diagonalize S:

        vec s;
        mat U;
        eig_sym(s, U, S);

        mat V = U*diagmat(1.0/sqrt(s));

        F = V.t() * F * V;


        vec eps;
        mat Cmat;
        eig_sym(eps, Cmat, F);

        C = V*Cmat.col(0);
        C = normalize(C);

        cout << C << endl;
        /*-----------------------------------------------------------------------------------------------------------*/
        //Calculate energy:

        double Eg=0.0;

        for(uint p=0; p < nOrbitals; p++){
            for(uint q=0; q < nOrbitals; q++){
                Eg += C(p)*C(q)*h(p,q);
            }
        }

        Eg = 2*Eg;
        for(uint p=0; p < nOrbitals; p++){
            for(uint r=0; r < nOrbitals; r++){
                for(uint q=0; q< nOrbitals; q++){
                    for(uint s=0; s < nOrbitals; s++){
                        Eg +=Q[p][r][q][s]*C(p)*C(q)*C(r)*C(s);
                    }
                }
            }
        }

        cout <<"Energy: " << Eg <<" step: " << step << endl;

    }




    /*-----------------------------------------------------------------------------------------------------------*/

    //De-Allocate memory to prevent memory leak
    for (uint i = 0; i < nOrbitals; ++i) {
        for (uint j = 0; j < nOrbitals; ++j){
            for (uint k = 0; k < nOrbitals; ++k){
                delete [] Q[i][j][k];
            }
            delete [] Q[i][j];
        }
        delete [] Q[i];
    }
    delete [] Q;

    return 0;
}


/*-----------------------------------------------------------------------------------------------------------*/
double nuclearRepulsion(){


    return 1;
}

/*-----------------------------------------------------------------------------------------------------------*/
vec normalize(vec C){

    double normFactor = norm(C,2);

    return C/normFactor;

}

/*-----------------------------------------------------------------------------------------------------------*/
double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   rowvec Rp,rowvec Rr,rowvec Rq,rowvec Rs,
                                   vec alpha, mat S){


    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];

    rowvec Ra = (alpha[p]*Rp + alpha[q]*Rq)/(alpha[p]+alpha[q]);
    rowvec Rb = (alpha[r]*Rr + alpha[s]*Rs)/(alpha[r]+alpha[s]);


    double t = (alpha[p]+alpha[q])*(alpha[r]+alpha[s])/ (alpha[p]+alpha[q]+alpha[r]+alpha[s]) * dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*errorFunction(t)* S(p,q)*S(s,r);

    return arg;


}


/*-----------------------------------------------------------------------------------------------------------*/
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

/*-----------------------------------------------------------------------------------------------------------*/
double nuclearAttarctionIntegral(const double p, const double q, const int Z, const int Rp, const int Rq, const mat R){

    double factor = 1.0/(p+q);
    double Rpq =dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-p*q*factor*Rpq);

    rowvec Rmc = (p*R.row(Rp) + q*R.row(Rq))*factor;


    double F0p = errorFunction(1.0/factor*dot(Rmc-R.row(0),Rmc-R.row(0)));
    double F0q = errorFunction(1.0/factor*dot(Rmc-R.row(1),Rmc-R.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0p+F0q);

    return nucAtt;

}

/*-----------------------------------------------------------------------------------------------------------*/
double kineticIntegral(const double p, const double q, rowvec Rp, rowvec Rq){

    double factor  = p*q/(p+q);
    double Rpq     = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-factor*Rpq);
    double kin     = 0.5*factor*(6-4*factor*Rpq)*pow(acos(-1)/(p+q),3.0/2.0)*expTerm;

    return kin;
}


/*-----------------------------------------------------------------------------------------------------------*/
double overlapIntegral(const double p, const double q, rowvec Rp, rowvec Rq){

    double factor = 1.0/(p+q);
    double Rpq = dot(Rp-Rq,Rp-Rq);
    double expTerm = exp(-p*q*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;

}


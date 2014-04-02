#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>

#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;


double overlapIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);
double overlapIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);

double kineticIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);
double kineticIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);

double nuclearAttarctionIntegral(const int Z, const uint p, const uint q,
                                 const uint Rp, const uint Rq, vec alpha, mat R);
double nuclearAttarctionIntegral_derivative(const int Z, const uint p, const uint q,
                                            const uint Rp, const uint Rq, vec alpha, mat R);

double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   const int  Rp, const int  Rr, const int  Rq, const int  Rs,
                                   vec alpha,mat R, mat S);
double electronInteractionIntegral_derivative(const int p, const int r, const int q, const int s,
                                              const int  Rp, const int  Rr, const int  Rq, const int  Rs,
                                              vec alpha, mat R);

double nuclearRepulsion(const mat R);
double nuclearRepulsion_derivative(const mat R);

double errorFunction(double arg);
double errorFunction_derivative(double arg);

double calculateEnergy(double ****Q, const mat h, const vec C, const mat R);
double calculateEnergy_derivative(double ****dQ, const mat dh, const vec C, const mat R);

vec normalize(vec C, mat S);

void writeToFile(const mat R, int n);

int main()
{
    boost::mpi::environment env;
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();
    /*-----------------------------------------------------------------------------------------------------------*/

    vector<Atom *> atoms;
    atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_0-4G.tm", { -0.675, 0, 0 }));
    atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_0-4G.tm", { 0.675, 0, 0 }));
    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);
    HFsolver* solver;
    solver = new RHF(system);
//    solver->setupOneParticleMatrix();
//    solver->setupTwoParticleMatrix();
    solver->runSolver();


    //System configuration:
    uint nBasisFunc = 4;
    uint nNuclei    = 2;
    int Z           = 1;

    uint e_nSteps  = 500;
    uint n_nSteps  = 1;

    double e_dt    = 0.1;
    double n_dt    = 4.3;

    double e_gamma = 1.0;
    double n_gamma = 0.0;

    double mu = 1.0;
    double M  = 1836.5*mu*0.5;

    double Eg,dEg;
    double a, b, c;
    double lambda;


    vec alpha = zeros(nBasisFunc);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;

    /*-----------------------------------------------------------------------------------------------------------*/
    //Initilize:
    uint nOrbitals = nBasisFunc * nNuclei;

    mat h = zeros(nOrbitals,nOrbitals);
    mat G = zeros(nOrbitals,nOrbitals);
    mat S = zeros(nOrbitals,nOrbitals);
    mat F = zeros(nOrbitals,nOrbitals);


    mat dh = zeros(nOrbitals,nOrbitals);
    mat dS = zeros(nOrbitals,nOrbitals);

    mat R  = zeros(nNuclei,3);

    double X, Xplus, Xminus;

    vec C      = ones(nOrbitals)*0.125;
    vec Cminus = ones(nOrbitals)*0.125;
    vec Cplus  = zeros(nOrbitals);

    vec2 lambdaVec;

    double ****Q;
    double ****dQ;
    Q = new double***[nOrbitals];
    dQ = new double***[nOrbitals];
    for (uint i = 0; i < nOrbitals; ++i) {
        Q[i]  = new double**[nOrbitals];
        dQ[i] = new double**[nOrbitals];

        for (uint j = 0; j < nOrbitals; ++j){
            Q[i][j]  = new double*[nOrbitals];
            dQ[i][j] = new double*[nOrbitals];

            for (uint k = 0; k < nOrbitals; ++k){
                Q[i][j][k]  = new double[nOrbitals];
                dQ[i][j][k] = new double[nOrbitals];
            }
        }
    }


    //One nuclei at -0.5ex and the other at 0.5ex
    X = 1.35;
    Xminus  = X;
    R(0,0)  = -X*0.5;
    R(1,0)  =  X*0.5;
    /*-----------------------------------------------------------------------------------------------------------*/

    for(uint nStep = 0; nStep < n_nSteps; nStep++){

        //Set up the h and S matrix:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rq = 0; Rq < R.n_rows; Rq++){

                for(uint p=0; p <alpha.size(); p++){
                    for(uint q=0; q <alpha.size(); q++){

                        S(p+Rp*4,q+Rq*4) = overlapIntegral(p,q,Rp,Rq,alpha,R);
                        h(p+Rp*4,q+Rq*4) = kineticIntegral(p,q,Rp,Rq,alpha,R)
                                + nuclearAttarctionIntegral(Z,p,q,Rp,Rq,alpha,R);
                    }
                }
            }

        }

//        cout << "S: "  << endl<< S
//             << "h: "  << endl<< h << endl;
//        sleep(5);
        C = normalize(C,S);
        Cminus = normalize(Cminus,S);
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

                                        Q[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral(p,r,q,s,Rp,Rr,Rq,Rs,alpha,R,S);
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }


        /*-----------------------------------------------------------------------------------------------------------*/

        //Loop over time
        for(uint eStep=0; eStep <= e_nSteps; eStep++){

            //Zero out elements
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
            F = h;

            const mat* F_ = solver->fockMatrix()[0];
//            cout << F << (*F_) << endl;
//            sleep(7);
            /*-----------------------------------------------------------------------------------------------------------*/
            //Calculate Ctilde:
            Cplus = (2*mu*C - (mu-e_gamma*e_dt*0.5)*Cminus - 4*e_dt*e_dt*F * C)/(mu+e_gamma*e_dt*0.5);


//            mat C_plus = (2*mu*C - (mu-e_gamma*e_dt*0.5)*C_minus - 4*e_dt*e_dt*F * C)/(mu+e_gamma*e_dt*0.5);


            //Calculate lambda:
            a = dot(C, S*S*S*C);
            b = -2*dot(S*Cplus,S*C);
            c = dot(Cplus,S*Cplus)-1.0;


            if(b*b - 4*a*c < 0 ){
                cerr << "Complex!" <<endl;
                cerr << a << "   "<< b << "   "<< c <<endl;
                exit (EXIT_FAILURE);

            }

            lambdaVec(0)  = (-b - sqrt(b*b - 4*a*c)) /(2*a);
            lambdaVec(1)  = (-b + sqrt(b*b - 4*a*c)) /(2*a);

            if(lambdaVec(1) < 0 && lambdaVec(0) < 0 ){
                cerr << "negative roots!!" <<endl;
                cerr << lambdaVec <<endl;
                exit (EXIT_FAILURE);

            }

            if(lambdaVec(0) >= 0.0 ){
                lambda = lambdaVec(0);
            }else{
                lambda = lambdaVec(1) ;
            }

            //Calculate C(t+h):
            Cplus -= lambda*S*C;


            //update C
            Cminus = C;
            C = Cplus;

            /*-----------------------------------------------------------------------------------------------------------*/
            //Calculate energy:
            Eg = calculateEnergy(Q,h,C,R);
//            cout << Eg << "," << endl;
            /*-----------------------------------------------------------------------------------------------------------*/


            //            cout.precision(8);
            //            cout <<"Energy: " << Eg <<" step: " << eStep << endl;
        }// Endof time loop electrons



        //Set up the dh and dS matrix:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rq = 0; Rq < R.n_rows; Rq++){

                for(uint p=0; p <alpha.size(); p++){
                    for(uint q=0; q <alpha.size(); q++){

                        dS(p+Rp*4,q+Rq*4) = overlapIntegral_derivative(p,q,Rp,Rq,alpha,R);
                        dh(p+Rp*4,q+Rq*4) = kineticIntegral_derivative(p,q,Rp,Rq,alpha,R)
                                + nuclearAttarctionIntegral_derivative(Z,p,q,Rp,Rq,alpha,R);
                    }
                }
            }

        }


        /*-----------------------------------------------------------------------------------------------------------*/
        //Set up the dQ array:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rr = 0; Rr < R.n_rows; Rr++){
                for(uint Rq = 0; Rq < R.n_rows; Rq++){
                    for(uint Rs = 0; Rs < R.n_rows; Rs++){

                        for(uint p=0; p <alpha.size(); p++){
                            for(uint r=0; r <alpha.size(); r++){
                                for(uint q=0; q <alpha.size(); q++){
                                    for(uint s=0; s <alpha.size(); s++){

                                        dQ[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral_derivative(p,r,q,s,Rp,Rr,Rq,Rs,alpha,R);
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

        /*-----------------------------------------------------------------------------------------------------------*/
        //Calculate dE/dX
        dEg = calculateEnergy_derivative(dQ,dh,C,R);

        /*-----------------------------------------------------------------------------------------------------------*/

        lambda =(mu+e_gamma*e_dt*0.5)/(2*e_dt*e_dt)*lambda;
        Xplus = (2*X*M - Xminus*(M-n_gamma*n_dt*0.5) - n_dt*n_dt*(dEg + lambda*dot(C,dS*C)) )/(M + n_gamma*n_dt*0.5 );

        cout.precision(8);
        cout << "[" << X << "," << Eg << "],"<< endl;
//            cout << X << "," << endl;

        R(0,0) = -Xplus*0.5;
        R(1,0) =  Xplus*0.5;

        writeToFile(R,nStep);

        Xminus = X;
        X      = Xplus;


    }// End of time loop nuclei
    /*-----------------------------------------------------------------------------------------------------------*/

    //De-Allocate memory to prevent memory leak
    for (uint i = 0; i < nOrbitals; ++i) {
        for (uint j = 0; j < nOrbitals; ++j){
            for (uint k = 0; k < nOrbitals; ++k){
                delete [] Q[i][j][k];
                delete [] dQ[i][j][k];
            }
            delete [] Q[i][j];
            delete [] dQ[i][j];
        }
        delete [] Q[i];
        delete [] dQ[i];
    }
    delete [] Q;
    delete [] dQ;

    return 0;
}
/*################################################################################################################*/
/*End of main function
##################################################################################################################*/


void writeToFile(const mat R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/data"<< n <<".xyz";
    myfile.open(outName.str().c_str(),ios::binary);
    myfile << 2    << "\n";
    myfile << "Hydrogen atoms  " << "\n";

    for(uint i=0;  i < R.n_rows; i++){
        myfile << R.row(i);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}

/*-----------------------------------------------------------------------------------------------------------*/
double calculateEnergy(double ****Q, const mat h, const vec C,const mat R){
    double Eg = 0.0;


    //One-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint q=0; q < C.n_elem; q++){
            Eg += C(p)*C(q)*h(p,q);
        }
    }
    Eg = 2*Eg;

    //Two-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint r=0; r < C.n_elem; r++){
            for(uint q=0; q< C.n_elem; q++){
                for(uint s=0; s < C.n_elem; s++){
                    Eg +=Q[p][r][q][s]*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }



    //Nuclear repulsion term
        Eg +=nuclearRepulsion(R);

    return Eg;

}
/*-----------------------------------------------------------------------------------------------------------*/
double calculateEnergy_derivative(double ****dQ, const mat dh, const vec C,const mat R){
    double dEg = 0.0;


    //One-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint q=0; q < C.n_elem; q++){
            dEg += C(p)*C(q)*dh(p,q);
        }
    }
    dEg = 2*dEg;

    //Two-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint r=0; r < C.n_elem; r++){
            for(uint q=0; q< C.n_elem; q++){
                for(uint s=0; s < C.n_elem; s++){
                    dEg +=dQ[p][r][q][s]*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }

    //    Nuclear repulsion term
    dEg +=nuclearRepulsion_derivative(R);

    return dEg;

}

/*-----------------------------------------------------------------------------------------------------------*/
vec normalize(vec C, mat S){
    double normFactor= 0.0;

    for(uint i= 0; i < C.n_elem; i++){
        for(uint j= 0; j < C.n_elem; j++){
            normFactor += C(i)*S(i,j)*C(j);
        }
    }
    return C/sqrt(normFactor);
}
/*-----------------------------------------------------------------------------------------------------------*/
double errorFunction(double t){

    if (t < 1.0E-6){
        return 1.0;
    }

    else{
        t = sqrt(t);
        double f = 1.0/t * erf(t)*sqrt(acos(-1))/2.0;
        return f;
    }

}

/*-----------------------------------------------------------------------------------------------------------*/
double overlapIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor = 1.0/(alpha(p)+ alpha(q));
    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-alpha(p)*alpha(q)*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;

}
/*-----------------------------------------------------------------------------------------------------------*/
double kineticIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-factor*Rpq);
    double kin     = 0.5*factor*(6-4*factor*Rpq)*pow(acos(-1)/(alpha(p)+ alpha(q)),3.0/2.0)*expTerm;

    return kin;
}
/*-----------------------------------------------------------------------------------------------------------*/
double nuclearAttarctionIntegral(const int Z, const uint p, const uint q,
                                 const uint Rp, const uint Rq, vec alpha, mat R){

    double factor = 1.0/(alpha(p)+ alpha(q));
    double Rpq =dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-alpha(p)*alpha(q)*factor*Rpq);

    rowvec Rmc = (alpha(p)*R.row(Rp) + alpha(q)*R.row(Rq))*factor;

    double F0p = errorFunction(1.0/factor*dot(Rmc-R.row(0),Rmc-R.row(0)));
    double F0q = errorFunction(1.0/factor*dot(Rmc-R.row(1),Rmc-R.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0p+F0q);

    return nucAtt;

}
/*-----------------------------------------------------------------------------------------------------------*/
double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   const int Rp, const int Rr, const int Rq, const int Rs,
                                   vec alpha, mat R, mat S){


    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];

    rowvec Ra = (alpha[p]*R.row(Rp) + alpha[q]*R.row(Rq))/A;
    rowvec Rb = (alpha[r]*R.row(Rr) + alpha[s]*R.row(Rs))/B;


    double t = (A*B/(A + B))*dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*errorFunction(t)*S(p+Rp*4,q+Rq*4)*S(s+Rs*4,r+Rr*4);

    return arg;


}
/*-----------------------------------------------------------------------------------------------------------*/
double nuclearRepulsion(const mat R){


    return 1/sqrt(dot(R.row(0) - R.row(1),R.row(0)- R.row(1)));
}

/*-----------------------------------------------------------------------------------------------------------*/
double nuclearRepulsion_derivative(const mat R){

    return -1/dot(R.row(0) - R.row(1),R.row(0)- R.row(1));
}





/*-----------------------------------------------------------------------------------------------------------*/
double errorFunction_derivative(double t){

    if (t < 1.0E-6){
        return -1.0/3.0;
    }

    else{
        double f = errorFunction(t);
        return (exp(-t)-f)/(2*t) ;
    }

}


//*-----------------------------------------------------------------------------------------------------------*/
double overlapIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));

    if(Rpq==0){
        return 0.0;
    }
    else{
        double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
        double X       = sqrt(dot(R.row(0)-R.row(1),R.row(0)-R.row(1)));
        double Spq     = overlapIntegral(p, q, Rp,Rq,alpha,R); // Use precalculated???
        double overlap = -2*factor*X*Spq;

        return overlap;
    }

}
/*-----------------------------------------------------------------------------------------------------------*/
double kineticIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));

    if(Rpq==0){
        return 0.0;
    }
    else{
        double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
        double X       = sqrt(dot(R.row(0)-R.row(1),R.row(0)-R.row(1)));
        double dSdX    = overlapIntegral_derivative(p, q, Rp, Rq,alpha,R);
        double Spq     = overlapIntegral(p, q, Rp, Rq,alpha,R); // Use precalculated???
        double kin     = -4*factor*factor*X*Spq + (3*factor - 2*factor*factor*X*X)*dSdX;

        return kin;
    }
}



/*-----------------------------------------------------------------------------------------------------------*/
double nuclearAttarctionIntegral_derivative(const int Z, const uint p, const uint q,
                                            const uint Rp, const uint Rq, vec alpha, mat R){

    double pq    = (alpha(p)+ alpha(q));
    double X     = sqrt(dot(R.row(0)-R.row(1),R.row(0)-R.row(1)));
    double Rpq   = sqrt(dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq)));
    double theta = 2*sqrt(pq/acos(-1));
    double Spq   = overlapIntegral(p, q, Rp,Rq,alpha,R);
    double nucAtt;


    if(Rpq == 0){
        double t = pq*X*X;
        nucAtt   = 2*Z*theta*Spq*X*pq*errorFunction_derivative(t);
    }
    else{
        double dSdX    = overlapIntegral_derivative(p, q, Rp, Rq, alpha, R);

        double F0p = errorFunction(X*X*alpha[p]*alpha[p]/pq);
        double F0q = errorFunction(X*X*alpha[q]*alpha[q]/pq);

        double dF0p = errorFunction_derivative(X*X*alpha[p]*alpha[p]/pq);
        double dF0q = errorFunction_derivative(X*X*alpha[q]*alpha[q]/pq);

        nucAtt = Z*theta*dSdX*(F0p + F0q) + 2*Z*(theta/pq)*Spq *(dF0p*alpha[p]*alpha[p] + dF0q*alpha[q]*alpha[q])*X;

    }


    return -nucAtt;
}

/*-----------------------------------------------------------------------------------------------------------*/
double electronInteractionIntegral_derivative(const int p, const int r, const int q, const int s,
                                              const int Rp, const int Rr, const int Rq, const int Rs,
                                              vec alpha, mat R){

    double X = sqrt(dot(R.row(0)-R.row(1),R.row(0)-R.row(1)));
    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];

    rowvec P = (alpha[p]*R.row(Rp)+ alpha[q]*R.row(Rq))/A;
    rowvec Q = (alpha[r]*R.row(Rr)+ alpha[s]*R.row(Rs))/B;

    double rho = 2*sqrt((A*B)/(acos(-1)*(A+B)));
    double t     = dot(P-Q,P-Q)*(A*B)/(A+B);


    double dSpqdX    = overlapIntegral_derivative(p, q, Rp, Rq, alpha, R);
    double dSrsdX    = overlapIntegral_derivative(r, s, Rr, Rs, alpha, R);

    double Spq = overlapIntegral(p, q, Rp, Rq, alpha, R);
    double Srs = overlapIntegral(r, s, Rr, Rs, alpha, R);


    double F0  = errorFunction(t);
    double dF0 = errorFunction_derivative(t);

    double arg = rho*dSpqdX*Srs*F0
            + rho*dSrsdX*Spq*F0
            + rho*Spq*Srs*dF0
            * (A*B)/(A+B) * 2*dot(P-Q,P-Q)/X;


    return arg;


}

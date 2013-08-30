//  Note that you need a main function with various include statemets in c++, like
//  #include <cmath> etc.
using namespace std;
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;


// compute the integral part of the Coulomb matrix element: <12|V|34>
double coulomb(const int, const int, const int, const int, const int, const int, const int, const int);

//compute the expectation value of the one-body operatot
double onebody(const int n, const int ml,double hbar, double omega);

//compute (-1)^k
int minusPower(const int);

// computes log(n!)
double LogFac (const int);

// Computes the first ratio in the Asinimovas expression
double LogRatio1(const int,const int,const int,const int);

// Computes the 2nd ratio in the Asinimovas expression
double LogRatio2(const int);

// computes first product of indices in the Anisimovas/Matulis expression
double Product1 (const int, const int, const int, const int, const int, const int, const int, const int);

// Computes the log of the 2nd product in the Asinimovas expression
double LogProduct2(const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int,const int);

// Computes the log of the 3rd product in the Asinimovas expression
double LogProduct3(const int,const int,const int,const int,const int,const int,const int,const int);

// The function lgamma() computes the logarithm of the gamma function of real argument x
double lgamma(double x);


// Intialise:
void initialise(int& num_particles, double& hbar, double& omega, int& max_iterations,int& max_energy);

class SingleParticleState {
public:
    int n;
    int ml;
    int ms;
    int index;
    double energy;
};

int main()
{
    // Variables:
    int num_particles,max_iterations,max_energy;
    double hbar,omega;

    // Quantum numbers:
    vec n = {0,1,2,4,5};
    vec ml = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    vec ms= {-1, 1};


    vector<SingleParticleState*> listOfStates;
    vec eigval;
    mat eigvec;

    // initialization
    initialise(num_particles,hbar,omega,max_iterations, max_energy);


    //initialization of states:
    for (uint i=0; i<n.size(); i++){
        for(uint j=0; j<ml.size(); j++){
            for(uint k=0; k<ms.size(); k++){

                double E_onebody=onebody(n(i),ml(j),hbar, omega);
                if( E_onebody <= max_energy) {
                    SingleParticleState* state = new SingleParticleState();
                    state->energy = E_onebody;
                    state->ml = ml(j);
                    state->n = n(i);
                    state->ms = ms(k);
                    listOfStates.push_back(state);
                }
            }
        }
    }

    //initilization of C and h_HF matrix
    mat C =eye(listOfStates.size(),listOfStates.size());
    mat h_HF(listOfStates.size(),listOfStates.size());


    //sorting of possible
    for(uint i=0;i<listOfStates.size();i++){
        for(uint j=0;j<i;j++){
            if((listOfStates.at(j)->energy == listOfStates.at(i)->energy && listOfStates.at(j)->ml > listOfStates.at(i)->ml) ||
                    listOfStates.at(j)->energy > listOfStates.at(i)->energy){

                SingleParticleState* temp;
                temp=listOfStates.at(i);
                listOfStates.at(i)=listOfStates.at(j);
                listOfStates.at(j)=temp;
            }
        }
    }

    //HF-algorithm:
    double E_temp=0.0;
    for (int i=0; i<max_iterations; i++){
        for (uint n1=0; n1<listOfStates.size(); n1++){
            for (uint n3=0; n3<listOfStates.size(); n3++){

                SingleParticleState* alpha=listOfStates.at(n1);
                SingleParticleState* gamma=listOfStates.at(n3);

                //Calculating the onebody energy:
                if (n1 ==n3){
                    h_HF(n1,n3)=onebody(alpha->n,alpha->ml,hbar, omega);}
                else {
                    h_HF(n1,n3)=0;}

                //Calculating the twobody energy:
                double E_coulumb=0.0;
                for(int p=0; p<num_particles; p++){
                    for(uint n2=0; n2<listOfStates.size();n2++){
                        for(uint n4=0;n4<listOfStates.size();n4++){

                            SingleParticleState* beta=listOfStates.at(n2);
                            SingleParticleState* delta=listOfStates.at(n4);

                            //Calculating the direct- and exchange term
                            double direct = coulomb(alpha->n, alpha->ml,beta->n,beta->ml,delta->n,delta->ml,gamma->n,gamma->ml);
                            double exchange = coulomb(alpha->n, alpha->ml,beta->n,beta->ml,gamma->n,gamma->ml,delta->n,delta->ml);

                            //Coefficients
                            double Cs = C(p,n2)*C(p,n4);

                            //Including spin degrees of freedom
                            //for <alpha,beta|gamma,delta>_AS
                            if(alpha->ms == gamma->ms && beta->ms == delta->ms) {
                                E_coulumb += Cs*(direct);
                            }
                            if(alpha->ms == delta->ms && beta->ms ==gamma->ms ) {
                                E_coulumb -= Cs*(exchange);
                            }
                        }}}
                h_HF (n1,n3)+= E_coulumb;
            }
        }

        //Diagonalization of h_HF matrix
        //The two/six first rows of eigvec are the two/six eigenvectors with lowest eigenvalue
        eig_sym(eigval, eigvec, h_HF);
        C=trans(eigvec);

        //Calculating energy with our new coffesients:

        //The one-body term:
        double E=0;
        for (int a=0; a<num_particles; a++){
            for (uint n1=0; n1<listOfStates.size(); n1++){
                for (uint n2=0; n2<listOfStates.size(); n2++){
                    SingleParticleState* alpha=listOfStates.at(n1);
                    if (n1 ==n2){
                        E+=C(a,n1)*C(a,n2)*onebody(alpha->n,alpha->ml,hbar,omega);}
                }
            }
        }

        //The two-body term:
        for (int a=0; a<num_particles; a++){
            for(int b=0; b<num_particles; b++){
                for (uint n1=0; n1<listOfStates.size(); n1++){
                    for (uint n2=0; n2<listOfStates.size(); n2++){
                        for (uint n3=0; n3<listOfStates.size(); n3++){
                            for (uint n4=0; n4<listOfStates.size(); n4++){

                                double E_coulumb=0.0;
                                SingleParticleState* alpha=listOfStates.at(n1);
                                SingleParticleState* beta=listOfStates.at(n2);
                                SingleParticleState* gamma=listOfStates.at(n3);
                                SingleParticleState* delta=listOfStates.at(n4);

                                double direct = coulomb(alpha->n, alpha->ml,beta->n,beta->ml,delta->n,delta->ml,gamma->n,gamma->ml);
                                double exchange = coulomb(alpha->n, alpha->ml,beta->n,beta->ml,gamma->n,gamma->ml,delta->n,delta->ml);

                                //Including spin degrees of freedom
                                //for <alpha,beta|gamma,delta>_AS
                                if(alpha->ms == gamma->ms && beta->ms == delta->ms) {
                                    E_coulumb += (direct);
                                }
                                if(alpha->ms == delta->ms && beta->ms ==gamma->ms ) {
                                    E_coulumb -=(exchange);
                                }
                                E += 0.5*C(a,n1)*C(b,n2)*C(a,n3)*C(b,n4)*E_coulumb;
                            }
                        }

                    }
                }
            }
        }

        cout << "Iteration number: "<< i+1 << ", Energy: " << setprecision(9) << E << endl;

        //compare the single-particle energies from the previous iteration
        //with those obtained from the new diagonalization.
        //If the total diﬀerence is smaller than a preﬁxed value, the iterative process is stopped
        if (abs(E_temp-E) < 10e-5){
            break;
        }
        E_temp=E;
    }

    cout << "Number of states: " << listOfStates.size() << endl;

    return 0;
}



    void initialise(int& num_particles, double& hbar, double& omega, int& max_iterations, int& max_energy)
    {

        num_particles=6;
        max_iterations=20;
        max_energy=3;
        hbar=1;
        omega=1;
        cout << "------Hartree-Fock -----" << endl;
        cout <<"-------------------------"<<endl;
        cout << "Number of particles: " << num_particles << endl;
        cout << "Number of shells: " << max_energy << endl;

    }

    double onebody(const int n, const int ml, double hbar, double omega){
        double e_onebody=0.0;

        e_onebody= hbar*omega*(2*n+abs(ml)+1);

        return e_onebody;
    }

    // compute one non-antisymmetric part of the Coulomb matrix element: <12|V|34> [see Asinimovas and Matulis (1997) for closed form  expression]
    // does NOT compute <12|V|34>as.
    // Be careful, this function doesn't account for spin.
    //
    // ex.: exchangeT= \delta(ms1,ms4) * \delta(ms2,ms3) * coulomb(n1,ml1,n2,ml2,n3,ml3,n4,ml4);
    //      directT  = \delta(ms1,ms3) * \delta(ms2,ms4) * coulomb(n1,ml1,n2,ml2,n4,ml4,n3,ml3);
    //
    double coulomb (const int n1, const int m1, const int n2, const int m2, const int n3, const int m3, const int n4, const int m4)
    {
        double coulombint=0.0;
        if(m1+m2 == m3+m4)
        {
            double temp;
            int lambda;

            int gamma1=0;
            int gamma2=0;
            int gamma3=0;
            int gamma4=0;
            int G=0;
            for(int j1=0;j1<=n1;j1++)
                for(int j2=0;j2<=n2;j2++)
                    for(int j3=0;j3<=n3;j3++)
                        for(int j4=0;j4<=n4;j4++)
                        {
                            gamma1 = (int) (j1 +j4 + 0.5*(abs(m1)+m1) + 0.5*(abs(m4)-m4));
                            gamma2 = (int) (j2 +j3 + 0.5*(abs(m2)+m2) + 0.5*(abs(m3)-m3));
                            gamma3 = (int) (j2 +j3 + 0.5*(abs(m3)+m3) + 0.5*(abs(m2)-m2));
                            gamma4 = (int) (j1 +j4 + 0.5*(abs(m4)+m4) + 0.5*(abs(m1)-m1));

                            G=gamma1+gamma2+gamma3+gamma4;

                            double Lgratio1 = LogRatio1(j1,j2,j3,j4);
                            double Lgproduct2 = LogProduct2(n1,m1,n2,m2,n3,m3,n4,m4,j1,j2,j3,j4);
                            double Lgratio2 = LogRatio2(G);

                            temp=0.0;
                            lambda=0;
                            for(int l1=0;l1<=gamma1;l1++)
                                for(int l2=0;l2<=gamma2;l2++)
                                    for(int l3=0;l3<=gamma3;l3++)
                                        for(int l4=0;l4<=gamma4;l4++)
                                        {
                                            lambda=l1+l2+l3+l4;
                                            if( (l1+l2)==(l3+l4) )
                                                temp += minusPower(gamma2+gamma3-l2-l3)*exp(LogProduct3(l1,l2,l3,l4,gamma1,gamma2,gamma3,gamma4)+lgamma(1+lambda*0.5)+lgamma((G-lambda+1)*0.5));
                                        }


                            coulombint += minusPower(j1+j2+j3+j4) *exp(Lgratio1 + Lgproduct2 + Lgratio2) *temp;

                        }
            coulombint *= Product1(n1,m1,n2,m2,n3,m3,n4,m4);
        }
        return coulombint;
    }


    // Computes the first ratio in the Asinimovas expression
    double LogRatio1(const int j1,const int j2,const int j3,const int j4)
    {
        double temp = -LogFac(j1)-LogFac(j2)-LogFac(j3)-LogFac(j4);
        return temp;
    }

    // Computes the 2nd ratio in the Asinimovas expression
    double LogRatio2(const int G)
    {
        double temp = -1*(G+1)*0.5*log(2);
        return temp;
    }


    // Computes the log of the 2nd product in the Asinimovas expression
    double LogProduct2(const int n1,const int m1,const int n2,const int m2,const int n3,const int m3,const int n4,const int m4,const int j1,const int j2,const int j3,const int j4)
    {
        double temp = LogFac(n1+abs(m1))+LogFac(n2+abs(m2))+LogFac(n3+abs(m3))+LogFac(n4+abs(m4))-LogFac(n1-j1)-LogFac(n2-j2)-LogFac(n3-j3)-LogFac(n4-j4)-LogFac(j1+abs(m1))-LogFac(j2+abs(m2))-LogFac(j3+abs(m3))-LogFac(j4+abs(m4));
        return temp;
    }


    // Computes the log of the 3rd product in the Asinimovas expression
    double LogProduct3(const int l1,const int l2,const int l3,const int l4,const int gamma1,const int gamma2,const int gamma3,const int gamma4)
    {
        double temp = LogFac(gamma1)+LogFac(gamma2)+LogFac(gamma3)+LogFac(gamma4)-LogFac(l1)-LogFac(l2)-LogFac(l3)-LogFac(l4)-LogFac(gamma1-l1)-LogFac(gamma2-l2)-LogFac(gamma3-l3)-LogFac(gamma4-l4);
        return temp;
    }


    // computes (-1)^k
    int minusPower(const int k)
    {
        int temp=abs(k%2)*-2+1;// gives 1 if k is even, gives -1 if k is odd
        return temp;
    }


    // computes log(n!)
    double LogFac (const int n)
    {
        if(n>170)
        {
            printf("#### Too big integer in LogFac(n)!!!!#### \n\n");
            exit(1);
        }
        double temp=0;
        for(int i=2; i<n+1;i++)
            temp+=log(i);
        return temp;
    }

    // computes first product of indices in the Anisimovas/Matulis expression
    double Product1 (const int n1, const int ml1, const int n2, const int ml2, const int n3, const int ml3, const int n4, const int ml4)
    {
        double temp=0;
        temp = LogFac(n1)+LogFac(n2)+LogFac(n3)+LogFac(n4)-LogFac(n1+abs(ml1))-LogFac(n2+abs(ml2))-LogFac(n3+abs(ml3))-LogFac(n4+abs(ml4));
        temp*=0.5;
        return exp(temp);
    }

    // lgamma.cpp -- log gamma function of real argument.
    //      Algorithms and coefficient values from "Computation of Special
    //      Functions", Zhang and Jin, John Wiley and Sons, 1996.
    //
    //  (C) 2003, C. Bond. All rights reserved.
    //
    //  Returns log(gamma) of real argument.
    //  NOTE: Returns 1e308 if argument is 0 or negative.
    //  taken on the web [http://www.crbond.com/math.htm]
    double lgamma(double x)
    {
#if CHECK
        if(x>=171)
        {
            printf("#### Too big integer in lgamma(x) to give accurate result!!!!#### \n\n");
            exit(1);
        }
#endif

        double x0,x2,xp,gl,gl0;
        int n,k;
        static double a[] = {
            8.333333333333333e-02,
            -2.777777777777778e-03,
            7.936507936507937e-04,
            -5.952380952380952e-04,
            8.417508417508418e-04,
            -1.917526917526918e-03,
            6.410256410256410e-03,
            -2.955065359477124e-02,
            1.796443723688307e-01,
            -1.39243221690590};

        x0 = x;
        if (x <= 0.0) return 1e308;
        else if ((x == 1.0) || (x == 2.0)) return 0.0;
        else if (x <= 7.0) {
            n = (int)(7-x);
            x0 = x+n;
        }
        x2 = 1.0/(x0*x0);
        xp = 2.0*M_PI;
        gl0 = a[9];
        for (k=8;k>=0;k--) {
            gl0 = gl0*x2 + a[k];
        }
        gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
        if (x <= 7.0) {
            for (k=1;k<=n;k++) {
                gl -= log(x0-1.0);
                x0 -= 1.0;
            }
        }
        return gl;
    }// end of lgamma function




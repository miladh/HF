#include "boys.h"

Boys::Boys()
{
}

double Boys::boysFunction(double arg, uint n)
{
    double F;

    if(arg <= 0.09){
        F = taylorExpandendBoys(arg,n,6);

    }else if (arg >= 20){
        F = doubleFactorial(2*n - 1) / pow(2,n+1) * sqrt(M_PI / pow(arg,2*n+1));
    }else{
        F = 0.0;
    }

    //        cout <<setprecision(10) << F <<endl;
    return F;


    //    if(arg < 0){
    //        throw("evaluating Boys function with negative argument!");

    //    }else if(arg == 0){
    //        F = 1.0 / (2*n + 1);

    //    }else{


}


double Boys::taylorExpandendBoys(double arg, uint n, int nterms){

    double F = 0.0;

    for(int k = 0; k < nterms; k++){
        F += pow(-arg,k) / (factorial(k) * (2*n + 2*k + 1));
    }

    return F;

}

double Boys::downwardRecursion(double arg, uint n){

    // p = previous, n = next
    double Fn = boysFunction(arg, n);
    double Fp = (2* arg * Fn + exp(-arg)) / (2*n -1);

    return Fp;

}

int Boys::doubleFactorial(int n)
{
    double result=1.0;
    for(int i=n; i>=1; i-=2){result *=i;}
    return result;
}

int Boys::factorial(int n)
{
    double result = 1;
    for (int i = 1; i <= n; i++){result *=i;}
    return result;
}



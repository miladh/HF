#include "boys.h"

Boys::Boys(uint highestOrder):
    m_highestOrder(highestOrder),
    m_results(zeros<rowvec>(highestOrder+1))
{
    readBoysForSmallArguments();
}

void Boys::readBoysForSmallArguments(){
    string path = "../../hf/infiles/tabulatedBoys.dat";
    m_Ftabulated.load(path);
}

void Boys::evaluateBoysFunctions(const double &arg)
{
    m_arg = arg;
    results();
}

rowvec Boys::getBoysFunctions() const
{
    return m_results;
}

void Boys::results()
{
    if(m_arg <= 50){
        m_results(m_highestOrder) = taylorExpandendBoys();
    }else{
        m_results(m_highestOrder) = asymptoticBoys();
    }

    downwardRecursion();
}


double Boys::taylorExpandendBoys(uint nterms, double dxt) const
{
    double FnMax = 0.0;
    int index = floor((m_arg + 0.5*dxt) / dxt);
    double dx = m_arg - index * dxt;

    for(uint k = 0; k < nterms; k++){
        FnMax  += (m_Ftabulated(index, m_highestOrder + k )*pow(-dx, k)) / factorial(k) ;
    }
    return FnMax;

}

double Boys::asymptoticBoys() const
{
    double FnMax;
    FnMax = doubleFactorial(2*m_highestOrder - 1) / pow(2,m_highestOrder+1)
            * sqrt(M_PI / pow(m_arg,2*m_highestOrder+1));
    return FnMax;
}


void Boys::downwardRecursion(){

    for(uint n = m_highestOrder; n > 0; n-- ){
        m_results(n-1) = (2* m_arg * m_results(n)  + exp(-m_arg)) / (2*n - 1);
    }
}

double Boys::doubleFactorial(const int &n)
{
    double result = 1.0;
    for(int i = n; i >= 1; i -= 2){result *=i;}
    return result;
}


double Boys::factorial(const int &n)
{
    double result = 1.0;
    for(int i = 1; i <= n; i++){result *=i;}
    return result;
}



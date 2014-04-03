#include "electronrepulsionintegral.h"

using namespace hf;
ElectronRepulsionIntegral::ElectronRepulsionIntegral(const int highestOrder,
                                                     const field<cube> *Eab,
                                                     const field<cube> *Ecd,
                                                     const PrimitiveGTO *primitiveA,
                                                     const PrimitiveGTO *primitiveB,
                                                     const PrimitiveGTO *primitiveC,
                                                     const PrimitiveGTO *primitiveD):

    m_R(new HermiteIntegrals(highestOrder)),
    m_Eab(Eab),
    m_Ecd(Ecd),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB),
    m_primitiveC(primitiveC),
    m_primitiveD(primitiveD)

{
//        cout << "--------------"<< endl;
//        cout << "Electron: "  << endl;
    //    cout << "coreC adr   "  << m_sourceCharge<< endl;
    //    cout << "primA adr   "  << m_primitiveB<< endl;
    //    cout << "primB adr   "  << m_primitiveA<< endl;

//        cout <<"E adr:   " << &m_Eab << endl;
//        cout <<"E1 adr:   " << (m_Eab->at(0))(0,0,0) <<"    " << &m_Eab->at(0) << endl;
//        cout <<"E2 adr:   " << &(m_Eab->at(1)) << endl;
//        cout <<"E3 adr:   " << &(m_Eab->at(2)) << endl;

//        cout << "--------------"<< endl;

//        cout <<"E adr:   " << &m_Ecd << endl;
//        cout <<"E1 adr:   " << (m_Ecd->at(0))(0,0,0) <<"    " << &m_Ecd->at(0) << endl;
//        cout <<"E2 adr:   " << &(m_Ecd->at(1)) << endl;
//        cout <<"E3 adr:   " << &(m_Ecd->at(2)) << endl;

}

double ElectronRepulsionIntegral::evaluate()
{
    const rowvec &A = m_primitiveA->center();
    const rowvec &B = m_primitiveB->center();
    const rowvec &C = m_primitiveC->center();
    const rowvec &D = m_primitiveD->center();

    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();
    const double &c  = m_primitiveC->exponent();
    const double &d  = m_primitiveD->exponent();

    double p = a + b;
    double q = c + d;

    double alpha = p*q/(p+q);
    rowvec PQ = (a*A + b*B)/p - (c*C + d*D)/q;


    double result = 0.0;
    int tMax = m_primitiveA->xPower() + m_primitiveB->xPower() + 1;
    int uMax = m_primitiveA->yPower() + m_primitiveB->yPower() + 1;
    int vMax = m_primitiveA->zPower() + m_primitiveB->zPower() + 1;
    int kMax = m_primitiveC->xPower() + m_primitiveD->xPower() + 1;
    int lMax = m_primitiveC->yPower() + m_primitiveD->yPower() + 1;
    int mMax = m_primitiveC->zPower() + m_primitiveD->zPower() + 1;

    m_R->updateR(PQ,alpha,
                 tMax + kMax - 2, uMax + lMax - 2, vMax + mMax - 2);


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){

                double Etuv = m_Eab->at(0)(m_primitiveA->xPower() , m_primitiveB->xPower(), t)
                            * m_Eab->at(1)(m_primitiveA->yPower() , m_primitiveB->yPower(), u)
                            * m_Eab->at(2)(m_primitiveA->zPower() , m_primitiveB->zPower(), v);

                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            double Eklm = m_Ecd->at(0)(m_primitiveC->xPower() , m_primitiveD->xPower(), k)
                                        * m_Ecd->at(1)(m_primitiveC->yPower() , m_primitiveD->yPower(), l)
                                        * m_Ecd->at(2)(m_primitiveC->zPower() , m_primitiveD->zPower(), m);
                            result += Etuv * Eklm * m_R->R(0, t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                        }
                    }
                }

            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q))
            * m_primitiveA->weight() * m_primitiveB->weight()
            * m_primitiveC->weight() * m_primitiveD->weight();

    return result;
}

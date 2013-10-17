#include "integrator.h"

Integrator::Integrator()
{
    rowvec RA = {1.2, 2.3, 3.4};
    R.row(0) = RA;
    rowvec RB = {-1.3, 1.4, -2.4};
    R.row(1) = RB;

    index(0,1)=1;
    primitives.push_back(new PrimitiveGTO(0.2, 1.0, "s1"));
    primitives.push_back(new PrimitiveGTO(0.3, 1.0, "s1"));

    tmax = l+3;

    for(int cor = 0; cor < 3; cor++){
        E[cor] = zeros(tmax, tmax, tmax);
    }

}

void Integrator::setE()
{
    double p = primitives[0]->alpha()+primitives[1]->alpha();
    double mu = primitives[0]->alpha()*primitives[1]->alpha()/p;

    rowvec P = (primitives[0]->alpha()*R.row(0) +  primitives[1]->alpha()*R.row(1))/p;
    rowvec Rab= R.row(0) - R.row(1);
    rowvec Rpa= P - R.row(0);
    rowvec Rpb= P - R.row(1);
    rowvec Kab = exp(-mu*Rab%Rab);


    for(int cor = 0; cor < 3; cor++){
        E[cor](0,0,0) = Kab(cor);
    }

    for(int cor=0; cor < 3; cor++){ //Loop for x,y,z

        for(int ib = 1; ib < tmax; ib++){
            for(int t = 0; t < tmax-1; t++){
                double E_0_ibp_tp = 0.0;
                double E_0_ibp_tn = 0.0;
                double E_0_ibp_tc = 0.0;

                if(!(t-1 < 0 || t-1 > (0 + (ib-1))) ) {

                    E_0_ibp_tp = E[cor](0,ib-1,t-1);
                }

                if(!(t+1 < 0 || t+1 > (0 + (ib-1))) ) {

                    E_0_ibp_tn = E[cor](0,ib-1,t+1);
                }

                if(!(t < 0 || t > (0 + (ib-1))) ) {

                    E_0_ibp_tc = E[cor](0,ib-1,t);
                }

                E[cor](0,ib,t) = 1.0/(2*p) * E_0_ibp_tp + Rpb(cor)*E_0_ibp_tc + (t+1)*E_0_ibp_tn;
            }
        }

        for(int ia = 1; ia < tmax; ia++){
            for(int ib = 0; ib < tmax; ib++){
                for(int t = 0; t < tmax-1; t++){
                    double E_iap_ibc_tp = 0.0;
                    double E_iap_ibc_tn = 0.0;
                    double E_iap_ibc_tc = 0.0;

                    if(!(t-1 < 0 || t-1 > (ia-1 + ib)) ) {

                        E_iap_ibc_tp = E[cor](ia-1,ib,t-1);
                    }

                    if(!(t+1 < 0 || t+1 > (ia-1 + ib)) ) {

                        E_iap_ibc_tn = E[cor](ia-1,ib,t+1);
                    }

                    if(!(t < 0 || t > (ia-1 + ib)) ) {

                        E_iap_ibc_tc = E[cor](ia-1,ib,t);
                    }

                    E[cor](ia,ib,t) = 1.0/(2*p) * E_iap_ibc_tp + Rpa(cor)*E_iap_ibc_tc + (t+1)*E_iap_ibc_tn;
                }
            }
        }



    }//End of cor=(x,y,z) loop


    cout <<"p " << p << endl;
    cout <<"mu "<< mu << endl;
    cout <<"P"<< P << endl;
    cout <<"ab"<< Rab <<endl;
    cout <<"pa"<< Rpa <<endl;
    cout <<"pb"<< Rpb <<endl;
    cout <<"kab"<< Kab <<endl;


    cout << E[0] << endl;
    cout << "-----------------------------------" <<endl;
    cout << E[1] << endl;
    cout << "-----------------------------------" <<endl;
    cout << E[2] << endl;




}

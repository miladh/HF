#ifndef DEFINES_H
#define DEFINES_H

#include <boost/regex.hpp>
#include <boost/config.hpp>

#define HFSOLVERTOLERANCE 1.0E-8
#define PROTONMASS 1836
#define time0 2.42E-17
#define E0 27.2114

#define USE_MPI 1
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)


class AtomMeta
{
public:
    int type;
    int charge;
    int nElectrons;
    double mass;


    static AtomMeta getData(std::string symbol){

        std::map <std::string, AtomMeta> atoms;
        AtomMeta atom;

        atom.type = 0;	   atom.charge = 0;	   atom.nElectrons = 0;	   atom.mass = 0;	             atoms["unknown"]= atom;
        atom.type = 1;	   atom.charge = 1;	   atom.nElectrons = 1;	   atom.mass = 1.0079;	        atoms["H"] =  atom;
        atom.type = 2;	   atom.charge = 2;	   atom.nElectrons = 2;	   atom.mass = 4.0026;	        atoms["He"] =  atom;
        atom.type = 3;	   atom.charge = 3;	   atom.nElectrons = 3;	   atom.mass = 6.941;	         atoms["Li"] =  atom;
        atom.type = 4;	   atom.charge = 4;	   atom.nElectrons = 4;	   atom.mass = 9.0122;	        atoms["Be"] =  atom;
        atom.type = 5;	   atom.charge = 5;	   atom.nElectrons = 5;	   atom.mass = 10.811;	        atoms["B"] =  atom;
        atom.type = 6;	   atom.charge = 6;	   atom.nElectrons = 6;	   atom.mass = 12.0107;	       atoms["C"] =  atom;
        atom.type = 7;	   atom.charge = 7;	   atom.nElectrons = 7;	   atom.mass = 14.0067;	       atoms["N"] =  atom;
        atom.type = 8;	   atom.charge = 8;	   atom.nElectrons = 8;	   atom.mass = 15.9994;	       atoms["O"] =  atom;
        atom.type = 9;	   atom.charge = 9;	   atom.nElectrons = 9;	   atom.mass = 18.9984;	       atoms["F"] =  atom;
        atom.type = 10;	   atom.charge = 10;	   atom.nElectrons = 10;	   atom.mass = 20.1797;	     atoms["Ne"] =  atom;
        atom.type = 11;	   atom.charge = 11;	   atom.nElectrons = 11;	   atom.mass = 22.9897;	     atoms["Na"] =  atom;
        atom.type = 12;	   atom.charge = 12;	   atom.nElectrons = 12;	   atom.mass = 24.305;	      atoms["Mg"] =  atom;
        atom.type = 13;	   atom.charge = 13;	   atom.nElectrons = 13;	   atom.mass = 26.9815;	     atoms["Al"] =  atom;
        atom.type = 14;	   atom.charge = 14;	   atom.nElectrons = 14;	   atom.mass = 28.0855;	     atoms["Si"] =  atom;
        atom.type = 15;	   atom.charge = 15;	   atom.nElectrons = 15;	   atom.mass = 30.9738;	     atoms["P"] =  atom;
        atom.type = 16;	   atom.charge = 16;	   atom.nElectrons = 16;	   atom.mass = 32.065;	      atoms["S"] =  atom;
        atom.type = 17;	   atom.charge = 17;	   atom.nElectrons = 17;	   atom.mass = 35.453;	      atoms["Cl"] =  atom;
        atom.type = 18;	   atom.charge = 18;	   atom.nElectrons = 18;	   atom.mass = 39.948;	      atoms["Ar"] =  atom;
        atom.type = 19;	   atom.charge = 19;	   atom.nElectrons = 19;	   atom.mass = 39.0983;	     atoms["K"] =  atom;
        atom.type = 20;	   atom.charge = 20;	   atom.nElectrons = 20;	   atom.mass = 40.078;	      atoms["Ca"] =  atom;
        atom.type = 21;	   atom.charge = 21;	   atom.nElectrons = 21;	   atom.mass = 44.9559;	     atoms["Sc"] =  atom;
        atom.type = 22;	   atom.charge = 22;	   atom.nElectrons = 22;	   atom.mass = 47.867;	      atoms["Ti"] =  atom;
        atom.type = 23;	   atom.charge = 23;	   atom.nElectrons = 23;	   atom.mass = 50.9415;	     atoms["V"] =  atom;
        atom.type = 24;	   atom.charge = 24;	   atom.nElectrons = 24;	   atom.mass = 51.9961;	     atoms["Cr"] =  atom;
        atom.type = 25;	   atom.charge = 25;	   atom.nElectrons = 25;	   atom.mass = 54.938;	      atoms["Mn"] =  atom;
        atom.type = 26;	   atom.charge = 26;	   atom.nElectrons = 26;	   atom.mass = 55.845;	      atoms["Fe"] =  atom;
        atom.type = 27;	   atom.charge = 27;	   atom.nElectrons = 27;	   atom.mass = 58.9332;	     atoms["Co"] =  atom;
        atom.type = 28;	   atom.charge = 28;	   atom.nElectrons = 28;	   atom.mass = 58.6934;	     atoms["Ni"] =  atom;
        atom.type = 29;	   atom.charge = 29;	   atom.nElectrons = 29;	   atom.mass = 63.546;	      atoms["Cu"] =  atom;
        atom.type = 30;	   atom.charge = 30;	   atom.nElectrons = 30;	   atom.mass = 65.39;	       atoms["Zn"] =  atom;
        atom.type = 31;	   atom.charge = 31;	   atom.nElectrons = 31;	   atom.mass = 69.723;	      atoms["Ga"] =  atom;
        atom.type = 32;	   atom.charge = 32;	   atom.nElectrons = 32;	   atom.mass = 72.64;	       atoms["Ge"] =  atom;
        atom.type = 33;	   atom.charge = 33;	   atom.nElectrons = 33;	   atom.mass = 74.9216;	     atoms["As"] =  atom;
        atom.type = 34;	   atom.charge = 34;	   atom.nElectrons = 34;	   atom.mass = 78.96;	       atoms["Se"] =  atom;
        atom.type = 35;	   atom.charge = 35;	   atom.nElectrons = 35;	   atom.mass = 79.904;	      atoms["Br"] =  atom;
        atom.type = 36;	   atom.charge = 36;	   atom.nElectrons = 36;	   atom.mass = 83.8;	        atoms["Kr"] =  atom;
        atom.type = 37;	   atom.charge = 37;	   atom.nElectrons = 37;	   atom.mass = 85.4678;	     atoms["Rb"] =  atom;
        atom.type = 38;	   atom.charge = 38;	   atom.nElectrons = 38;	   atom.mass = 87.62;	       atoms["Sr"] =  atom;
        atom.type = 39;	   atom.charge = 39;	   atom.nElectrons = 39;	   atom.mass = 88.9059;	     atoms["Y"] =  atom;
        atom.type = 40;	   atom.charge = 40;	   atom.nElectrons = 40;	   atom.mass = 91.224;	      atoms["Zr"] =  atom;
        atom.type = 41;	   atom.charge = 41;	   atom.nElectrons = 41;	   atom.mass = 92.9064;	     atoms["Nb"] =  atom;
        atom.type = 42;	   atom.charge = 42;	   atom.nElectrons = 42;	   atom.mass = 95.94;	       atoms["Mo"] =  atom;
        atom.type = 43;	   atom.charge = 43;	   atom.nElectrons = 43;	   atom.mass = 98;	          atoms["Tc"] =  atom;
        atom.type = 44;	   atom.charge = 44;	   atom.nElectrons = 44;	   atom.mass = 101.07;	      atoms["Ru"] =  atom;
        atom.type = 45;	   atom.charge = 45;	   atom.nElectrons = 45;	   atom.mass = 102.9055;	    atoms["Rh"] =  atom;
        atom.type = 46;	   atom.charge = 46;	   atom.nElectrons = 46;	   atom.mass = 106.42;	      atoms["Pd"] =  atom;
        atom.type = 47;	   atom.charge = 47;	   atom.nElectrons = 47;	   atom.mass = 107.8682;	    atoms["Ag"] =  atom;
        atom.type = 48;	   atom.charge = 48;	   atom.nElectrons = 48;	   atom.mass = 112.411;	     atoms["Cd"] =  atom;
        atom.type = 49;	   atom.charge = 49;	   atom.nElectrons = 49;	   atom.mass = 114.818;	     atoms["In"] =  atom;
        atom.type = 50;	   atom.charge = 50;	   atom.nElectrons = 50;	   atom.mass = 118.71;	      atoms["Sn"] =  atom;
        atom.type = 51;	   atom.charge = 51;	   atom.nElectrons = 51;	   atom.mass = 121.76;	      atoms["Sb"] =  atom;
        atom.type = 52;	   atom.charge = 52;	   atom.nElectrons = 52;	   atom.mass = 127.6;	       atoms["Te"] =  atom;
        atom.type = 53;	   atom.charge = 53;	   atom.nElectrons = 53;	   atom.mass = 126.9045;	    atoms["I"] =  atom;
        atom.type = 54;	   atom.charge = 54;	   atom.nElectrons = 54;	   atom.mass = 131.293;	     atoms["Xe"] =  atom;
        atom.type = 55;	   atom.charge = 55;	   atom.nElectrons = 55;	   atom.mass = 132.9055;	    atoms["Cs"] =  atom;
        atom.type = 56;	   atom.charge = 56;	   atom.nElectrons = 56;	   atom.mass = 137.327;	     atoms["Ba"] =  atom;
        atom.type = 57;	   atom.charge = 57;	   atom.nElectrons = 57;	   atom.mass = 138.9055;	    atoms["La"] =  atom;
        atom.type = 58;	   atom.charge = 58;	   atom.nElectrons = 58;	   atom.mass = 140.116;	     atoms["Ce"] =  atom;
        atom.type = 59;	   atom.charge = 59;	   atom.nElectrons = 59;	   atom.mass = 140.9077;	    atoms["Pr"] =  atom;
        atom.type = 60;	   atom.charge = 60;	   atom.nElectrons = 60;	   atom.mass = 144.24;	      atoms["Nd"] =  atom;
        atom.type = 61;	   atom.charge = 61;	   atom.nElectrons = 61;	   atom.mass = 145;	         atoms["Pm"] =  atom;
        atom.type = 62;	   atom.charge = 62;	   atom.nElectrons = 62;	   atom.mass = 150.36;	      atoms["Sm"] =  atom;
        atom.type = 63;	   atom.charge = 63;	   atom.nElectrons = 63;	   atom.mass = 151.964;	     atoms["Eu"] =  atom;
        atom.type = 64;	   atom.charge = 64;	   atom.nElectrons = 64;	   atom.mass = 157.25;	      atoms["Gd"] =  atom;
        atom.type = 65;	   atom.charge = 65;	   atom.nElectrons = 65;	   atom.mass = 158.9253;	    atoms["Tb"] =  atom;
        atom.type = 66;	   atom.charge = 66;	   atom.nElectrons = 66;	   atom.mass = 162.5;	       atoms["Dy"] =  atom;
        atom.type = 67;	   atom.charge = 67;	   atom.nElectrons = 67;	   atom.mass = 164.9303;	    atoms["Ho"] =  atom;
        atom.type = 68;	   atom.charge = 68;	   atom.nElectrons = 68;	   atom.mass = 167.259;	     atoms["Er"] =  atom;
        atom.type = 69;	   atom.charge = 69;	   atom.nElectrons = 69;	   atom.mass = 168.9342;	    atoms["Tm"] =  atom;
        atom.type = 70;	   atom.charge = 70;	   atom.nElectrons = 70;	   atom.mass = 173.04;	      atoms["Yb"] =  atom;
        atom.type = 71;	   atom.charge = 71;	   atom.nElectrons = 71;	   atom.mass = 174.967;	     atoms["Lu"] =  atom;
        atom.type = 72;	   atom.charge = 72;	   atom.nElectrons = 72;	   atom.mass = 178.49;	      atoms["Hf"] =  atom;
        atom.type = 73;	   atom.charge = 73;	   atom.nElectrons = 73;	   atom.mass = 180.9479;	    atoms["Ta"] =  atom;
        atom.type = 74;	   atom.charge = 74;	   atom.nElectrons = 74;	   atom.mass = 183.84;	      atoms["W"] =  atom;
        atom.type = 75;	   atom.charge = 75;	   atom.nElectrons = 75;	   atom.mass = 186.207;	     atoms["Re"] =  atom;
        atom.type = 76;	   atom.charge = 76;	   atom.nElectrons = 76;	   atom.mass = 190.23;	      atoms["Os"] =  atom;
        atom.type = 77;	   atom.charge = 77;	   atom.nElectrons = 77;	   atom.mass = 192.217;	     atoms["Ir"] =  atom;
        atom.type = 78;	   atom.charge = 78;	   atom.nElectrons = 78;	   atom.mass = 195.078;	     atoms["Pt"] =  atom;
        atom.type = 79;	   atom.charge = 79;	   atom.nElectrons = 79;	   atom.mass = 196.9665;	    atoms["Au"] =  atom;
        atom.type = 80;	   atom.charge = 80;	   atom.nElectrons = 80;	   atom.mass = 200.59;	      atoms["Hg"] =  atom;
        atom.type = 81;	   atom.charge = 81;	   atom.nElectrons = 81;	   atom.mass = 204.3833;	    atoms["Tl"] =  atom;
        atom.type = 82;	   atom.charge = 82;	   atom.nElectrons = 82;	   atom.mass = 207.2;	       atoms["Pb"] =  atom;
        atom.type = 83;	   atom.charge = 83;	   atom.nElectrons = 83;	   atom.mass = 208.9804;	    atoms["Bi"] =  atom;
        atom.type = 84;	   atom.charge = 84;	   atom.nElectrons = 84;	   atom.mass = 209;	         atoms["Po"] =  atom;
        atom.type = 85;	   atom.charge = 85;	   atom.nElectrons = 85;	   atom.mass = 210;	         atoms["At"] =  atom;
        atom.type = 86;	   atom.charge = 86;	   atom.nElectrons = 86;	   atom.mass = 222;	         atoms["Rn"] =  atom;
        atom.type = 87;	   atom.charge = 87;	   atom.nElectrons = 87;	   atom.mass = 223;	         atoms["Fr"] =  atom;
        atom.type = 88;	   atom.charge = 88;	   atom.nElectrons = 88;	   atom.mass = 226;	         atoms["Ra"] =  atom;
        atom.type = 89;	   atom.charge = 89;	   atom.nElectrons = 89;	   atom.mass = 227;	         atoms["Ac"] =  atom;
        atom.type = 90;	   atom.charge = 90;	   atom.nElectrons = 90;	   atom.mass = 232.0381;	    atoms["Th"] =  atom;
        atom.type = 91;	   atom.charge = 91;	   atom.nElectrons = 91;	   atom.mass = 231.0359;	    atoms["Pa"] =  atom;
        atom.type = 92;	   atom.charge = 92;	   atom.nElectrons = 92;	   atom.mass = 238.0289;	    atoms["U"] =  atom;
        atom.type = 93;	   atom.charge = 93;	   atom.nElectrons = 93;	   atom.mass = 237;	         atoms["Np"] =  atom;
        atom.type = 94;	   atom.charge = 94;	   atom.nElectrons = 94;	   atom.mass = 244;	         atoms["Pu"] =  atom;
        atom.type = 95;	   atom.charge = 95;	   atom.nElectrons = 95;	   atom.mass = 243;	         atoms["Am"] =  atom;
        atom.type = 96;	   atom.charge = 96;	   atom.nElectrons = 96;	   atom.mass = 247;	         atoms["Cm"] =  atom;
        atom.type = 97;	   atom.charge = 97;	   atom.nElectrons = 97;	   atom.mass = 247;	         atoms["Bk"] =  atom;
        atom.type = 98;	   atom.charge = 98;	   atom.nElectrons = 98;	   atom.mass = 251;	         atoms["Cf"] =  atom;
        atom.type = 99;	   atom.charge = 99;	   atom.nElectrons = 99;	   atom.mass = 252;	         atoms["Es"] =  atom;
        atom.type = 100;    atom.charge = 100;	   atom.nElectrons = 100;	   atom.mass = 257;	       atoms["Fm"] =  atom;
        atom.type = 101;	   atom.charge = 101;	   atom.nElectrons = 101;	   atom.mass = 258;	       atoms["Md"] =  atom;
        atom.type = 102;	   atom.charge = 102;	   atom.nElectrons = 102;	   atom.mass = 259;	       atoms["No"] =  atom;
        atom.type = 103;	   atom.charge = 103;	   atom.nElectrons = 103;	   atom.mass = 262;	       atoms["Lr"] =  atom;
        atom.type = 104;	   atom.charge = 104;	   atom.nElectrons = 104;	   atom.mass = 261;	       atoms["Rf"] =  atom;
        atom.type = 105;	   atom.charge = 105;	   atom.nElectrons = 105;	   atom.mass = 262;	       atoms["Db"] =  atom;
        atom.type = 106;	   atom.charge = 106;	   atom.nElectrons = 106;	   atom.mass = 266;	       atoms["Sg"] =  atom;
        atom.type = 107;	   atom.charge = 107;	   atom.nElectrons = 107;	   atom.mass = 264;	       atoms["Bh"] =  atom;
        atom.type = 108;	   atom.charge = 108;	   atom.nElectrons = 108;	   atom.mass = 277;	       atoms["Hs"] =  atom;
        atom.type = 109;	   atom.charge = 109;	   atom.nElectrons = 109;	   atom.mass = 268;	       atoms["Mt"] =  atom;

        atom.type = 1;	   atom.charge = 1;	   atom.nElectrons = 0;	   atom.mass = 1.0079;	        atoms["H+"] =  atom;
        atom.type = 6;	   atom.charge = 6;	   atom.nElectrons = 7;	   atom.mass = 12.0107;	       atoms["C-"] =  atom;




        if(atoms.find(symbol) != atoms.end())
        {
            return atoms[symbol];
        }else{
            return atoms["unknown"];
        }
    }
};


#endif // DEFINES_H

#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <fstream>
#include <iostream>
#include <boost/mpi.hpp>
#include <hf.h>


using namespace std;
using namespace arma;
using namespace hf;

SUITE(GD){
//    TEST(HF)
//    {
//        /*
//     * test case:   geometrical derivative of energy
//     * system:      H2
//     * basis:       6-31G**
//     * method:      RHF
//     * source:
//     *      numerical differentiation of energy
//     * */

//        int myRank = 0;
//#if USE_MPI
//        boost::mpi::communicator world;
//        myRank = world.rank();
//#endif

//        if(myRank == 0){
//            cout << "system:    " << "H2" << endl;
//            cout << "method:    " << "RHF" << endl;
//            cout << "basis:     " << "6-31G**" << endl;
//        }

//        //Initializing the system
//        vector<Atom *> atoms;
//        atoms.push_back( new Atom("infiles/turbomole/atom_6_basis_6-31Gds.tm", {-0.5, 0.0, 0.0}));
//        atoms.push_back( new Atom("infiles/turbomole/atom_6_basis_6-31Gds.tm", { 0.5, 0.0, 0.0}));

//        ElectronicSystem *system = new ElectronicSystem ();
//        system->addAtoms(atoms);

//        RHF *solver = new RHF(system);
//        Atom* atomA = atoms.at(0);
//        Atom* atomB = atoms.at(1);


//        GeometricalDerivative* m_GD = new GeometricalDerivative(system, solver);

//        //Domain
//        vec bondLength = linspace(2.0, 2.2, 1);
//        vec energies = 0*bondLength;
//        vec gradients = 0*bondLength;



//        for(uint x = 0; x < bondLength.n_elem; x++){
//            rowvec X = {bondLength(x) , 0 ,0 };

//            atomA->setCorePosition(X * -0.5);
//            atomB->setCorePosition(X *  0.5);
//            solver->runSolver();
//            energies[x] = solver->energy();

//            mat gradE = m_GD->energyGradient();
//            gradients[x] = gradE(1,0);
//        }

//        if(myRank == 0){
//        for (int i = 0 ; i < signed(bondLength.n_elem); i++){
//            cout << "[" << bondLength(i)<< "  ,  "  << energies(i) << "  ,  " << gradients(i) << "]," << endl;
//        }
//    }

//    }



//    TEST(N2)
//    {
//        /*
//         * test case:   geometrical derivative of energy
//         * system:      N2
//         * basis:       STO-3G
//         * method:      RHF
//         * source:
//         *      numerical differentiation of energy at minimum
//         *      Modern Qunatum Chemisty, Szabo
//         *
//         * */

//        int myRank = 0;
//#if USE_MPI
//        boost::mpi::communicator world;
//        myRank = world.rank();
//#endif

//        if(myRank == 0){
//            cout << "system:    " << "N2" << endl;
//            cout << "method:    " << "RHF" << endl;
//            cout << "basis:     " << "STO-3G" << endl;
//        }

//        double bondLength = 2.143;

//        //Initializing the system
//        vector<Atom *> atoms;
//        atoms.push_back( new Atom("infiles/turbomole/atom_7_basis_STO-3G.tm", {0.0, 0.0, 0.0}));
//        atoms.push_back( new Atom("infiles/turbomole/atom_7_basis_STO-3G.tm", {0.0, 0.0, 0.0}));

//        ElectronicSystem *system = new ElectronicSystem ();
//        system->addAtoms(atoms);

//        RHF *solver = new RHF(system);
//        GeometricalDerivative* m_GD = new GeometricalDerivative(system, solver);

//        Atom* atomA = atoms.at(0);
//        Atom* atomB = atoms.at(1);

//        //Calculate gradient analytically
//        rowvec X = {bondLength , 0 ,0 };
//        atomA->setCorePosition(X * -0.5);
//        atomB->setCorePosition(X * 0.5);
//        solver->runSolver();
//        mat gradE = m_GD->energyGradient();
//        double gradient = gradE(1,0);

//        //Calculate gradient numerically
//        double h = 1.0E-8;
//        rowvec dx = {h , 0 ,0 };


//        atomA->setCorePosition((X-dx) * -0.5);
//        atomB->setCorePosition((X-dx) *  0.5);
//        solver->runSolver();
//        double Ep = solver->energy();

//        atomA->setCorePosition((X+dx) * -0.5);
//        atomB->setCorePosition((X+dx) *  0.5);
//        solver->runSolver();
//        double En = solver->energy();

//        atomA->setCorePosition((X-2.0*dx) * -0.5);
//        atomB->setCorePosition((X-2.0*dx) *  0.5);
//        solver->runSolver();
//        double Epp = solver->energy();

//        atomA->setCorePosition((X+2.0*dx) * -0.5);
//        atomB->setCorePosition((X+2.0*dx) *  0.5);
//        solver->runSolver();
//        double Enn = solver->energy();

////                double numericalGradient = (-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
//        double numericalGradient = (En - Ep) / (2.0 * h);
//        CHECK_CLOSE(numericalGradient, gradient, 1e-8);
//        cout << bondLength << "    " <<numericalGradient << "    " << gradient << endl;

//    }


        TEST(H2_RHF)
        {
            /*
             * test case:   geometrical derivative of energy
             * system:      H2
             * basis:       STO-3G
             * method:      RHF
             * source:
             *      numerical differentiation of energy
             * */

            int myRank = 0;
    #if USE_MPI
            boost::mpi::communicator world;
            myRank = world.rank();
    #endif

            if(myRank == 0){
                cout << "system:    " << "H2" << endl;
                cout << "method:    " << "RHF" << endl;
                cout << "basis:     " << "STO-3G" << endl;
            }

            //Initializing the system
            vector<Atom *> atoms;
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.5, 0.0, 0.0}));
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", { 0.5, 0.0, 0.0}));

            ElectronicSystem *system = new ElectronicSystem ();
            system->addAtoms(atoms);

            RHF *solver = new RHF(system);
            GeometricalDerivative* m_GD = new GeometricalDerivative(system, solver);

            Atom* atomA = atoms.at(0);
            Atom* atomB = atoms.at(1);


            //Domain
            vec bondLength = linspace(1.0, 4.5, 10);
            vec gradient   = 0*bondLength;
            vec numericalGradient = 0*bondLength;

            //Calculate gradient analytically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };

                atomA->setCorePosition(X * -0.5);
                atomB->setCorePosition(X * 0.5);
                solver->runSolver();
                mat gradE = m_GD->energyGradient();
                gradient(x) = gradE(1,0);
            }

            double h = 1.0E-7;
            //Calculate gradient numerically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };
                rowvec dx = {h , 0 ,0 };


                atomA->setCorePosition((X-dx) * -0.5);
                atomB->setCorePosition((X-dx) *  0.5);
                solver->runSolver();
                double Ep = solver->energy();

                atomA->setCorePosition((X+dx) * -0.5);
                atomB->setCorePosition((X+dx) *  0.5);
                solver->runSolver();
                double En = solver->energy();

                atomA->setCorePosition((X-2.0*dx) * -0.5);
                atomB->setCorePosition((X-2.0*dx) *  0.5);
                solver->runSolver();
                double Epp = solver->energy();

                atomA->setCorePosition((X+2.0*dx) * -0.5);
                atomB->setCorePosition((X+2.0*dx) *  0.5);
                solver->runSolver();
                double Enn = solver->energy();

                numericalGradient(x) = (-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
            }


            for(uint i = 0; i < numericalGradient.n_elem; i++){
                CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-8);

            }

        }


        TEST(H2_UHF)
        {
            /*
             * test case:   geometrical derivative of energy
             * system:      H2
             * basis:       STO-3G
             * method:      UHF
             * source:
             *      numerical differentiation of energy
             * */

            int myRank = 0;
    #if USE_MPI
            boost::mpi::communicator world;
            myRank = world.rank();
    #endif

            if(myRank == 0){
                cout << "system:    " << "H2" << endl;
                cout << "method:    " << "UHF" << endl;
                cout << "basis:     " << "STO-3G" << endl;
            }

            //Initializing the system
            vector<Atom *> atoms;
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-0.5, 0.0, 0.0}));
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", { 0.5, 0.0, 0.0}));

            ElectronicSystem *system = new ElectronicSystem ();
            system->addAtoms(atoms);

            UHF *solver = new UHF(system);
            solver->setDampingFactor(0.0);
            GeometricalDerivative* m_GD = new GeometricalDerivative(system, solver);

            Atom* atomA = atoms.at(0);
            Atom* atomB = atoms.at(1);

            //Domain
            vec bondLength = linspace(1.0, 4.5, 10);
            vec gradient   = 0*bondLength;
            vec numericalGradient = 0*bondLength;

            //Calculate gradient analytically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };

                atomA->setCorePosition(X * -0.5);
                atomB->setCorePosition(X * 0.5);
                solver->runSolver();
                mat gradE = m_GD->energyGradient();
                gradient(x) = gradE(1,0);
            }

            double h = 1.0E-7;
            //Calculate gradient numerically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };
                rowvec dx = {h , 0 ,0 };


                atomA->setCorePosition((X-dx) * -0.5);
                atomB->setCorePosition((X-dx) *  0.5);
                solver->runSolver();
                double Ep = solver->energy();

                atomA->setCorePosition((X+dx) * -0.5);
                atomB->setCorePosition((X+dx) *  0.5);
                solver->runSolver();
                double En = solver->energy();

                atomA->setCorePosition((X-2.0*dx) * -0.5);
                atomB->setCorePosition((X-2.0*dx) *  0.5);
                solver->runSolver();
                double Epp = solver->energy();

                atomA->setCorePosition((X+2.0*dx) * -0.5);
                atomB->setCorePosition((X+2.0*dx) *  0.5);
                solver->runSolver();
                double Enn = solver->energy();

                numericalGradient(x) = (-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
            }


            for(uint i = 0; i < numericalGradient.n_elem; i++){
                CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-8);
            }

        }


        TEST(H2_RHF_6_31Gds)
        {
            /*
             * test case:   geometrical derivative of energy
             * system:      H2
             * basis:       6-31G**
             * method:      RHF
             * source:
             *      numerical differentiation of energy
             * */

            int myRank = 0;
    #if USE_MPI
            boost::mpi::communicator world;
            myRank = world.rank();
    #endif

            if(myRank == 0){
                cout << "system:    " << "H2" << endl;
                cout << "method:    " << "RHF" << endl;
                cout << "basis:     " << "6-31G**" << endl;
            }

            //Initializing the system
            vector<Atom *> atoms;
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {-0.5, 0.0, 0.0}));
            atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { 0.5, 0.0, 0.0}));

            ElectronicSystem *system = new ElectronicSystem ();
            system->addAtoms(atoms);

            RHF *solver = new RHF(system);

            GeometricalDerivative* m_GD = new GeometricalDerivative(system, solver);

            Atom* atomA = atoms.at(0);
            Atom* atomB = atoms.at(1);

            //Domain
            vec bondLength = linspace(1.0, 4.5, 10);
            vec gradient   = 0*bondLength;
            vec numericalGradient = 0*bondLength;

            //Calculate gradient analytically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };

                atomA->setCorePosition(X * -0.5);
                atomB->setCorePosition(X * 0.5);
                solver->runSolver();
                mat gradE = m_GD->energyGradient();
                gradient(x) = gradE(1,0);
            }

            double h = 1.0E-7;
            //Calculate gradient numerically
            for(uint x = 0; x < bondLength.n_elem; x++){
                rowvec X = {bondLength(x) , 0 ,0 };
                rowvec dx = {h , 0 ,0 };


                atomA->setCorePosition((X-dx) * -0.5);
                atomB->setCorePosition((X-dx) *  0.5);
                solver->runSolver();
                double Ep = solver->energy();

                atomA->setCorePosition((X+dx) * -0.5);
                atomB->setCorePosition((X+dx) *  0.5);
                solver->runSolver();
                double En = solver->energy();

                atomA->setCorePosition((X-2.0*dx) * -0.5);
                atomB->setCorePosition((X-2.0*dx) *  0.5);
                solver->runSolver();
                double Epp = solver->energy();

                atomA->setCorePosition((X+2.0*dx) * -0.5);
                atomB->setCorePosition((X+2.0*dx) *  0.5);
                solver->runSolver();
                double Enn = solver->energy();

                numericalGradient(x) = (-Enn + 8. * En - 8. * Ep + Epp) / (12.0 * h);
            }


            for(uint i = 0; i < numericalGradient.n_elem; i++){
                CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-8);

            }

        }



}

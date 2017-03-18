[![Project Status: Inactive - The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive) 

The project has reached a stable, usable state but is no longer being actively developed.

HF
====
An implementation of the Hartree-Fock method for atoms and molecules, based on McMurchie-Davidson scheme using Gaussian basis functions. The input basis files read by the code are in the Turbomole format. These files are taken from the website [Basis Set Exchange](https://bse.pnl.gov/bse/portal)


To compile the following libraries are needed:
- Armadillo
- Libconfig
- HDF5
- Boost

Editing of the code is best done using [QtCreator](http://qt-project.org/downloads).

Example of config file:
```
#-------------------------------
# Example of configuration file
# for the HF program
#-------------------------------
chemicalSystem =
{
    name = "H2O";
    atoms =
    (
        {
        basis = "atom_8_basis_STO-3G.tm";
        position = [0.0, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [1.797, 0.0, 0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [-0.448, -1.740, 0.0];
        }
    )
    
    # Optional (if both are set to zero, there will be equal number of each spin):
    nSpinUpElectrons = 0;
    nSpinDownElectrons = 0;
}

solverSettings =
{
    # 0 = "Restricted HF"
    # 1 = "Unrestricted HF"
    method = 0;

    maxNumOfIteration = 10000;
    dampingFactor     = 0.5;


    DIISprocedureSettings =
    {
        # 0 = "off"
        # 1 = "on"
        useDIISprocedure = 0;

        iterationLimit = 20;
        nTerms         = 3;

    };


};

analysisSettings=
{
    # 0 = "off"
    # 1 = "on"
    
    saveResults = 1;
    outputFilePath = "/"

    saveEnergies = 1;
    dipoleMoment = 1;
    atomicPartialCharge = 1;
    chargeDensity = 0;
    electrostaticPotential = 0;

    densitySettings=
    {
        nGridPoints  = 200;
        minGridPoint = -8.;
        maxGridPoint =  8.;
    };

};

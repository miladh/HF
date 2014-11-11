HF
====
An implementation of the Hartree-Fock method.
Editing of the code is best done using QtCreator.

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
        velocity = [0.        ,  0.        ,  0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [1.797, 0.0, 0.0];
        velocity = [0.        ,  0.        ,  0.0];
        }
        ,
        {
        basis = "atom_1_basis_STO-3G.tm";
        position = [-0.448, -1.740, 0.0];
        velocity = [0.        ,  0.        ,  0.0];
        }
    )
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
    saveResults = 1;
    outputFilePath = "/home/milad/kurs/qmd/"

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

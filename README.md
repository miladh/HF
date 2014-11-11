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
    name = "CH5";

    atoms =
    (
        {
        basis = "atom_1_anion_basis_6-31ppGds.tm";
        position = [0.        ,  0.        ,  -2.79679446];
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
    method =0;

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
    outputFilePath = "/"

    saveEnergies = 1;
    dipoleMoment = 1;
    atomicPartialCharge = 1;
    chargeDensity = 1;
    electrostaticPotential = 0;


    densitySettings=
    {
        nGridPoints  = 200;
        minGridPoint = -8.;
        maxGridPoint =  8.;
    };

};

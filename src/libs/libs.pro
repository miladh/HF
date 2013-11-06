include(../../defaults.pri) 
TEMPLATE = lib 
TARGET = hartree-fock

SOURCES += \
    primitiveGTO/primitiveGTO.cpp \
    contractedGTO/contractedGTO.cpp \
    system/system.cpp \
    integrator/integrator.cpp \
    math/boys.cpp \
    math/hermitecoefficients.cpp \
    hfSolver/hfsolver.cpp \
    basisSet/basisset.cpp \
    basisSet/h_321g.cpp \
    basisSet/h_quadzeta.cpp \
    basisSet/o_321g.cpp \
    basisSet/h_sto6.cpp

HEADERS += \
    primitiveGTO/primitiveGTO.h \
    contractedGTO/contractedGTO.h \
    system/system.h \
    integrator/integrator.h \
    math/boys.h \
    math/hermitecoefficients.h \
    hfSolver/hfsolver.h \
    basisSet/basisset.h \
    basisSet/h_321g.h \
    basisSet/h_quadzeta.h \
    basisSet/o_321g.h \
    basisSet/h_sto6.h

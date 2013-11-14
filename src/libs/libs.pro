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
    basisSet/splitValence/h_321g.cpp \
    basisSet/h_quadzeta.cpp \
    basisSet/splitValence/o_321g.cpp \
    math/hermiteintegrals.cpp \
    basisSet/splitValence/li_321g.cpp \
    basisSet/splitValence/splitvalence.cpp \
    basisSet/splitValence/h_431g.cpp \
    basisSet/splitValence/o_431g.cpp \
    cpmd/cpmd.cpp

HEADERS += \
    primitiveGTO/primitiveGTO.h \
    contractedGTO/contractedGTO.h \
    system/system.h \
    integrator/integrator.h \
    math/boys.h \
    math/hermitecoefficients.h \
    hfSolver/hfsolver.h \
    basisSet/basisset.h \
    basisSet/splitValence/h_321g.h \
    basisSet/h_quadzeta.h \
    basisSet/splitValence/o_321g.h \
    math/hermiteintegrals.h \
    basisSet/splitValence/li_321g.h \
    basisSet/splitValence/splitvalence.h \
    basisSet/splitValence/h_431g.h \
    basisSet/splitValence/o_431g.h \
    cpmd/cpmd.h

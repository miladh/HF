include(../defaults.pri)


TEMPLATE = lib
TARGET = ../lib/hartree-fock


SOURCES += \
    primitiveGTO/primitiveGTO.cpp \
    contractedGTO/contractedGTO.cpp \
    system/system.cpp \
    integrator/integrator.cpp \
    math/boys.cpp \
    math/hermitecoefficients.cpp \
    hfSolver/hfsolver.cpp \
    basisSet/basisset.cpp \
    math/hermiteintegrals.cpp \
    cpmd/cpmd.cpp \
    bomd/bomd.cpp \
    geometricalDerivative/geometricalderivative.cpp \
    hfSolver/uhf.cpp \
    hfSolver/rhf.cpp \
    atom/atom.cpp \
    parser/turbomoleparser.cpp

HEADERS += \
    primitiveGTO/primitiveGTO.h \
    contractedGTO/contractedGTO.h \
    system/system.h \
    integrator/integrator.h \
    math/boys.h \
    math/hermitecoefficients.h \
    hfSolver/hfsolver.h \
    basisSet/basisset.h \
    math/hermiteintegrals.h \
    cpmd/cpmd.h \
    bomd/bomd.h \
    defines.h \
    geometricalDerivative/geometricalderivative.h \
    hfSolver/uhf.h \
    hfSolver/rhf.h \
    atom/atom.h \
    parser/turbomoleparser.h

OTHER_FILES += ../include/hf.h ../install/include/hf.h

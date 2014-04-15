include(../defaults.pri)


TEMPLATE = lib
TARGET = ../lib/hartree-fock


SOURCES += \
    primitiveGTO/primitiveGTO.cpp \
    contractedGTO/contractedGTO.cpp \
    integrator/integrator.cpp \
    math/boys.cpp \
    math/hermitecoefficients.cpp \
    hfSolver/hfsolver.cpp \
    math/hermiteintegrals.cpp \
    bomd/bomd.cpp \
    geometricalDerivative/geometricalderivative.cpp \
    hfSolver/uhf.cpp \
    hfSolver/rhf.cpp \
    atom/atom.cpp \
    parser/turbomoleparser.cpp \
    system/electronicsystem.cpp \
    analyser/analyser.cpp \
    integrator/overlap/overlapintegral.cpp \
    integrator/kinetic/kineticintegral.cpp \
    integrator/nuclearAttraction/nuclearattractionintegral.cpp \
    integrator/electronRepulsion/electronrepulsionintegral.cpp \
    integrator/overlap/overlapintegralgd.cpp \
    integrator/kinetic/kineticintegralgd.cpp \
    integrator/nuclearAttraction/nuclearattractionintegralgd.cpp \
    integrator/electronRepulsion/electronrepulsionintegralgd.cpp \
    integrator/dipole/dipoleintegral.cpp \
    outputManager/outputmanager.cpp

HEADERS += \
    primitiveGTO/primitiveGTO.h \
    contractedGTO/contractedGTO.h \
    integrator/integrator.h \
    math/boys.h \
    math/hermitecoefficients.h \
    hfSolver/hfsolver.h \
    math/hermiteintegrals.h \
    bomd/bomd.h \
    defines.h \
    geometricalDerivative/geometricalderivative.h \
    hfSolver/uhf.h \
    hfSolver/rhf.h \
    atom/atom.h \
    parser/turbomoleparser.h \
    system/electronicsystem.h \
    analyser/analyser.h \
    integrator/overlap/overlapintegral.h \
    integrator/kinetic/kineticintegral.h \
    integrator/nuclearAttraction/nuclearattractionintegral.h \
    integrator/electronRepulsion/electronrepulsionintegral.h \
    integrator/overlap/overlapintegralgd.h \
    integrator/kinetic/kineticintegralgd.h \
    integrator/nuclearAttraction/nuclearattractionintegralgd.h \
    integrator/electronRepulsion/electronrepulsionintegralgd.h \
    integrator/dipole/dipoleintegral.h \
    outputManager/outputmanager.h


OTHER_FILES += ../include/hf.h ../install/include/hf.h


!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}

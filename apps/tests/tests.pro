include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = hf_tests

SOURCES += main.cpp \
    hfsolverTests/solverTests.cpp \
    BoysFunctionTests/boysTests.cpp \
    gradientTests/gradientTests.cpp \
    integratorTests/overlapIntegral.cpp \
    integratorTests/kineticIntegral.cpp \
    integratorTests/coulombIntegral_1.cpp


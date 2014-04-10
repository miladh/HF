include(../../defaults.pri)
include(../apps_defaults.pri)

TARGET = hf_tests

SOURCES += main.cpp \
    BoysFunctionTests/boysTests.cpp \
    gradientTests/gradientTests.cpp \
    integratorTests/overlapIntegral.cpp \
    integratorTests/kineticIntegral.cpp \
    hfsolverTests/rhfTest.cpp \
    hfsolverTests/uhfTest.cpp \
    integratorTests/electronRepulsionIntegral.cpp \
    integratorTests/nuclearAttractionIntegral.cpp


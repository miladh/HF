TEMPLATE = app

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lunittest++

include(../defaults.pri)

SOURCES += main.cpp \
    ../src/integrator/integrator.cpp\
    ../src/math/hermitecoefficients.cpp \
    ../src/primitiveGTO/primitiveGTO.cpp \
    ../src/math/boys.cpp


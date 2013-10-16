TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../defaults.pri)

SOURCES += main.cpp \
    primitiveGTO/primitiveGTO.cpp \
    contractedGTO/contractedGTO.cpp \
    system/system.cpp \
    integrator/integrator.cpp

HEADERS += \
    primitiveGTO/primitiveGTO.h \
    contractedGTO/contractedGTO.h \
    system/system.h \
    integrator/integrator.h


include(../defaults.pri)

TARGET = hartree-fock-tests

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lunittest++ -L$$TOP_OUT_PWD/src/libs -lhartree-fock

SOURCES += main.cpp

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += src/main.cpp

HEADERS +=

LIBS += -llapack -larmadillo

QMAKE_CXXFLAGS += -std=c++0x

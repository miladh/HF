TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += src/main.cpp \
    src/contractedGTO/contractedGTO.cpp \
    src/primitiveGTO/primitiveGTO.cpp

HEADERS += \
    src/contractedGTO/contractedGTO .h \
    src/primitiveGTO/primitiveGTO.h

LIBS += -llapack -larmadillo

QMAKE_CXXFLAGS += -std=c++0x

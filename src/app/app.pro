include(../../defaults.pri) 
TEMPLATE = app 
SOURCES = main.cpp

LIBS += -L$$TOP_PWD/src/libs  -lhartree-fock

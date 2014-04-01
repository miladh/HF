TEMPLATE = app

LIBS += -lhdf5 -lhdf5_cpp
LIBS += -lunittest++ -L$$TOP_OUT_PWD/lib -lhartree-fock
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)

INCLUDEPATH += $$TOP_PWD/include

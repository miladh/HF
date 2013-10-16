INCLUDEPATH += $$PWD
LIBS += -llapack -larmadillo -lconfig++

QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x


release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}


# MPI Settings
QMAKE_CXX = ccache mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

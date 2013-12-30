LIBS += -llapack -larmadillo -lconfig++ -lboost_regex

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS


CURRENT_COMPILER = $$QMAKE_CXX
QMAKE_CXX = ccache $$CURRENT_COMPILER

# Directories
INCLUDEPATH += $$TOP_PWD/src/libs
SRC_DIR = $$TOP_PWD


#Copy infiles
copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

unix {
    LIBS += -lblas -llapack -larmadillo
} else {
LIBS += -L$$PWD/../../../../../Armadillo/ -larmadillo
LIBS += -L$$PWD/../../../../../Armadillo/examples/lib_win64 -llapack_win64_MT -lblas_win64_MT
INCLUDEPATH += $$PWD/../../../../../Armadillo/include
}

SOURCES += main.cpp

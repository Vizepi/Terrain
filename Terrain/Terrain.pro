#-------------------------------------------------
#
# Project created by QtCreator 2016-11-20T09:20:33
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
greaterThan(QT_MAJOR_VERSION, 4): CONFIG+= c++11

TARGET = Terrain
TEMPLATE = app


SOURCES += main.cpp\
    Terrain.cpp \
    TerrainBuilder.inl \
    Perlin.cpp \
    Matrix.inl

HEADERS  += \
    Terrain.h \
    TerrainBuilder.h \
    Vector.h \
    AABB.h \
    Perlin.h \
    Matrix.h \
    PI.h \
    Tree.h

FORMS    += MainWindow.ui

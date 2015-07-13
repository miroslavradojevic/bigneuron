
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = /home/miroslav/vaa3d/v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/common_lib/include

HEADERS	+= NeurAuto_plugin.h \
    toolbox.h \
    model2d.h \
    model3d.h
SOURCES	+= NeurAuto_plugin.cpp \
    toolbox.cpp \
    model2d.cpp \
    model3d.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp
SOURCES += $$VAA3DPATH/v3d_main/basic_c_fun/basic_4dimage_create.cpp

TARGET	= $$qtLibraryTarget(NeurAuto)
DESTDIR	= $$VAA3DPATH/bin/plugins/bigneuronhackathon/NeurAuto/

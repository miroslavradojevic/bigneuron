
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = /home/miroslav/vaa3d/v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/common_lib/include

HEADERS	+= NeuronFloater_plugin.h \
    nf_dialog.h \
    toolbox.h \
    model.h \
    node.h \
    tracer.h
SOURCES	+= NeuronFloater_plugin.cpp \
    toolbox.cpp \
    model.cpp \
    node.cpp \
    tracer.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp

TARGET	= $$qtLibraryTarget(NeuronFloater)
DESTDIR	= $$VAA3DPATH/bin/plugins/bigneuronhackathon/NeuronFloater/

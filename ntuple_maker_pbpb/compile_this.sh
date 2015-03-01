rootcint -f linkdefs.cxx -c `root-config --cflags` defs/linkdefs.h
g++ -m64 -Wall linkdefs.cxx -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --libs --glibs` class_def/JetAna.C class_def/HLT.C class_def/HiTree.C class_def/Tracks.C class_def/Skim.C  class_def/GenParticles.C copy_ntuple.C -o copy_ntuple.out


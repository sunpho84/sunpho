#!/bin/bash

#Setta il livello di ottimizzazione
OTTIM=-O4

#Attiva o meno il campo magnetico
#CAMPO_EST="-D mf"

#Attiva o meno il potenziale chimico di isospin reale
POTRE=" -D ficp" 

#Attiva le misure chirali
MISU=" -D chiral_meas" 

gfortran -x f77-cpp-input -o main main.f $OTTIM $CAMPO_EST $POTRE $MISU -Wall #-D isotropic -Ddebug1

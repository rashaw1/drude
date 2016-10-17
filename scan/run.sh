#!/bin/bash
for f in ./scan2/*.input
do
    ../drudeh "$f" ../basissets/standard.basis
done

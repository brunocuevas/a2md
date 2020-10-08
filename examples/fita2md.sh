#!/bin/bash

workon amd

echo "---"

cp "${A2MDTEST}/benzene/benzene.mol2" .
cp "${A2MDTEST}/benzene/benzene.g09.out" .
cp "${A2MDTEST}/benzene/benzene.wfn" .

cdens -i benzene.wfn -u .benzene.s1.csv -r 5 -p 0.25
cdens -i benzene.wfn -u .benzene.s2.csv -r 5 -p 0.10
cdens -i benzene.wfn -u .benzene.s3.csv -r 5 -p 0.01

python3 "${A2MD}/a2mdutils.py" update-mol2 --charges=NPA benzene.mol2 benzene.wfn benzene.g09.out benzene.n.mol2

cat .benzene.s?.csv | sed "s/,/  /g" | sed "/^x/d" > benzene.density.csv
rm .benzene.s?.csv

python3 "${A2MD}/a2mdrun.py" fit --cluster=rbf --verbose=2 benzene.n.mol2 benzene.density.csv
python3 "${A2MD}/a2mdrun.py" write-dx --output=benzene.dx --expand=3.0 --res=0.2 benzene.n.mol2 benzene.n.ppp

echo "---"

#!/bin/bash

workon amd

echo "---"

cp "${A2MDTEST}/benzene/benzene.mol2" .
echo "preparing a dft calculation for a methane molecule"
python3 "${A2MD}/a2mdutils.py" prepare-qm --method=wB97X --basis="6-31G(d)" --program=g09 benzene.mol2
mv benzene.g09.input benzene.dft.g09.input
echo "preparing a dft calculation + NPA analysis for a methane molecule"
python3 "${A2MD}/a2mdutils.py" prepare-qm --method=MP2 --population=NPA --basis="6-31G(d)" --program=g09 benzene.mol2
echo "preparing a dft calculation in orca"
python3 "${A2MD}/a2mdutils.py" prepare-qm --method=wB97X --basis="6-31G(d)" --program=orca benzene.mol2

echo "---"

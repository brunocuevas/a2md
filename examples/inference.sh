#!/bin/bash

workon amd

echo "---"

cp "${A2MDTEST}/benzene/benzene.mol2" .
python3 "${A2MD}/conformational.py" optimize --output=benzene.mol2 benzene.mol2
python3 "${A2MD}/a2mdpredict.py" predict-model --predictor=a2mdc --device=cpu --output=benzene.ml.ppp benzene.mol2
python3 "${A2MD}/a2mdrun.py" write-dx --output=benzene.ml.dx --expand=3.0 --res=0.2 benzene.mol2 benzene.ml.ppp

echo "---"

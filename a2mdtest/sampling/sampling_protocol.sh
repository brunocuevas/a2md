Cdens -i ${1}.wfn -p 1.500 -r 25 -u ${1}_s01.csv
Cdens -i ${1}.wfn -p 0.250 -r 25 -u ${1}_s02.csv
Cdens -i ${1}.wfn -p 0.150 -r 25 -u ${1}_s03.csv
Cdens -i ${1}.wfn -p 0.100 -r 25 -u ${1}_s04.csv
Cdens -i ${1}.wfn -p 0.010 -r 25 -u ${1}_s05.csv
Cdens -i ${1}.wfn -p 0.001 -r 25 -u ${1}_s06.csv
Cdens -i ${1}.wfn -d ${1}_mp2.dx -r 0.15 -s 3.0

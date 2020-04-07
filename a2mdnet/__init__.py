from pathlib import Path
LIBRARY_PATH = Path(__file__)

ELEMENT2NN = {1: 0, 6: 1, 7: 2, 8: 3}
SYMBOL2NN = {"H":0, "C":1, "N":2, "O":3}
ELEMENT2SYMBOL = {1:'H', 6:'C', 7:'N', 8:'O'}
ALLOWED_SPECIES = [1, 6, 7, 8]
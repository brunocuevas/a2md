from pathlib import Path
LIBRARY_PATH = Path(__file__).parent
AEV_PARAMETERS = LIBRARY_PATH / "params/aev_params.params"
ELEMENT2NN = {1: 0, 6: 1, 7: 2, 8: 3}
SYMBOL2NN = {"H":0, "C":1, "N":2, "O":3}
ELEMENT2SYMBOL = {1:'H', 6:'C', 7:'N', 8:'O'}
ALLOWED_SPECIES = [1, 6, 7, 8]

FUNCTION_NAMES2POSITION = {
    "CR":("core", None),
    "VR":("iso", 1),
    "CVR":("iso", 0),
    "hVR":("iso", 0),
    "B01":("aniso", 0),
    "B02":("aniso", 1)
}

MODELS = dict(
    a2mdc = LIBRARY_PATH / "models/a2mdnet_coeff.pt"
)
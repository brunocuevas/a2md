import torch
import json
import pathlib
from a2mdnet.modules import SymFeats
from a2mdnet.data import PolymerDataset
from a2mdnet import LIBRARY_PATH

NEW_AEVS = LIBRARY_PATH / "params/aev_polymer.params"

if __name__ == '__main__':

    print("--- sym feats test")
    print("--- testing S symmetry features")

    with open('F:/edip/dynamic/ccc/.ccc.input') as f:
        idx = json.load(f)[:10]

    pd = PolymerDataset(
        device=torch.device('cuda:0'), dtype=torch.float,
        ids=idx,
        molecular_data_path=pathlib.Path('F:/edip/dynamic/ccc/'),
        model_parameters_path=pathlib.Path('F:/edip/dynamic/ccc/'),
        prefetch=True
    )

    l, t, x, q, s, sq, ii, ai, ip, ap = pd[0]
    l.unsqueeze_(0)
    x.unsqueeze_(0)
    sf = SymFeats(NEW_AEVS)
    sf.to(torch.device('cuda:0'))
    u = sf.forward(l, x)
    print("DONE!")


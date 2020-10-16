import torch
from torch import nn
from a2mdio.params import AMDParameters
from a2mdnet.density_models import GenAMD
from a2mdnet.modules import QPChargeNormalization
import numpy as np


class AMDMolecule(nn.Module):

    def __init__(
            self, coordinates: np.ndarray, labels: np.ndarray, parameters: AMDParameters,
            charge: np.ndarray,
            device: torch.device, dtype: torch.dtype
    ):
        """

        :param coordinates:
        :param labels:
        :param parameters:
        :param device:
        :param dtype:
        """
        super(AMDMolecule, self).__init__()
        self.labels = torch.tensor(data=labels, dtype=torch.long, device=device)
        self.coordinates = torch.tensor(data=coordinates, dtype=dtype, device=device)
        self.model = GenAMD(parameters, device, dtype)
        self.normalizer = QPChargeNormalization(device=device)
        self.charge = torch.tensor(charge, dtype=dtype, device=device)
        n = coordinates.shape[0]
        m = parameters.get_maxfunctions()
        self.coefficients = nn.Parameter(
            data=0.1*torch.randn(n, m, dtype=dtype, device=device), requires_grad=True
        )

    def forward(self, sample):

        coefficients = self.coefficients.unsqueeze(0)
        coordinates = self.coordinates.unsqueeze(0)
        labels = self.labels.unsqueeze(0)
        charge = self.charge.unsqueeze(0)
        integrals = self.model.integrate(labels)
        norm_coefficients = self.normalizer.forward_iso(
            coefficients=coefficients, integrals=integrals, charges=charge
        )

        p = self.model.forward(sample, norm_coefficients, labels, coordinates)
        return p

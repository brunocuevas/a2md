from a2md.utils import integrate_from_dict, integrate_from_old_dict
from a2md.support import SupportAngular
import math

if __name__ == "__main__":

    u = integrate_from_dict(
        dict(
            bond= None,
            center= 0,
            coefficient= 0.545682446307695,
            params= {
                "A": 1/math.pi,
                "B": 2.0
            },
            support_type= "VR"
        )
    )

    v = integrate_from_old_dict(
        {
            "center": 19,
            "support_type": "AG",
            "coefficient": 0.19373057071284497,
            "params": {
                "Alpha": 3.5,
                "G": 2.5,
                "U": 1.29,
                "Psi":0
            },
            "bond": 21
        }
    )
    ax = SupportAngular(
        coordinates=[0.0, 0.0, 0.0],
        A=1.29, B=2.5, alpha=3.5, P=1
    )
    print(u)
    print(v)
    print(ax.integral())
    print("DONE!")
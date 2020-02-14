from a2mdlib.qm import CubeFile
import matplotlib.pyplot as plt
import numpy as np
if __name__ == "__main__":
    cf = CubeFile(file=r'C:/Users/Bruno/ownCloud/main/a2md/output/ep/MP2_lb/water.cube')
    cf.read()
    x = cf.get_plane(0, axis='x')

    plt.contourf(np.log(x))
    plt.show()
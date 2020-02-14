import numpy as np
from a2md.utils import Wavefunction


class sphere :
    def __init__(self, center, radius):
        """

        :param center:
        :type center: np.ndarray
        :param radius:
        :type radius: float
        """
        self.__center = center
        self.__radius = radius
    def __set_bins(self, phi_res, theta_res):
        """

        :param phi_res:
        :param theta_res:
        :return:
        """
        n_theta = int(self.__radius * theta_res)
        d_theta = np.pi / n_theta
        d_phis  = 2 * np.pi / phi_res
        bins = np.ones(n_theta +1, dtype='int64')
        for i in range(n_theta + 1) :
            per = 2 * np.pi * self.__radius * np.sin(i*d_theta)
            if not np.isclose(per, 0.0) :
                bins[i] = int(per / d_phis) + 1
        return n_theta, bins

    def estimate_points(self, phi_res, theta_res):
        n_points_theta, bins_phi = self.__set_bins(phi_res, theta_res)
        n_points_phi = np.sum(bins_phi) - 1
        return n_points_phi

    def trace(self, phi_res, theta_res):
        """
        Using the wide-adopted convention:
        - theta for the angle with z axis
        - phi for the angle with x axis projected on xy plane
        :param phi_res:
        :type phi_res: int
        :param theta_res:
        :type theta_res: int
        :return:
        """
        n_points_theta, res_by_parallel = self.__set_bins(phi_res, theta_res)
        total_points = np.sum(res_by_parallel)

        cartesian = np.zeros((total_points - 1, 3), dtype='float64')
        k = 0
        for i in range(n_points_theta):
            theta = i * np.pi / n_points_theta
            for j in range(res_by_parallel[i]):
                phi = (2 * np.pi) * float(j) / res_by_parallel[i]
                cartesian[k, 0] = np.cos(phi) * np.sin(theta)
                cartesian[k, 1] = np.sin(phi) * np.sin(theta)
                cartesian[k, 2] = np.cos(theta)
                k += 1
        return (cartesian * self.__radius) + self.__center


    def restrict(self, coords):
        """

        :param coords:
        :type coords: np ndarray
        :return:
        """
        d = np.linalg.norm(coords - self.__center, axis=1) + 1e-3
        T = d >= self.__radius

        return T

class surface :
    def __init__(self, center_coordinates, radius_list, factors = None):
        """

        :param center_coordinates:
        :param radius_list:
        """
        sphere_list = []
        natoms = center_coordinates.shape[0]
        if natoms != len(radius_list):
            raise IOError("length of list of atoms must be equal to the coordinates defined")
        if factors is None :
            factors = np.ones(len(sphere_list), dtype='float64')
        for center in range(natoms) :
            sphere_list.append(
                sphere(
                    center=center_coordinates[center,:],
                    radius=radius_list[center]*factors[center] / 25.0
                )
            )
        self.__sphere_list = sphere_list

    def calculate_surface_points(self, phi_resolution = 10, theta_resolution = 5):
        """
        Resolution is provided as sampled points by au of contour
        :param phi_resolution:
        :param theta_resolution:
        :return:
        """
        initial_surface_size = [
            sph.estimate_points(phi_resolution, theta_resolution) for sph in self.__sphere_list
        ]
        initial_surface = np.zeros((np.sum(initial_surface_size), 3), dtype='float64')
        prev = 0
        for i, sph in enumerate(self.__sphere_list) :
            post = initial_surface_size[i] + prev
            initial_surface[prev:post] = sph.trace(phi_resolution, theta_resolution)
            prev = prev + initial_surface_size[i]
        restriction = np.ones(np.sum(initial_surface_size), dtype='bool')
        for i, sph in enumerate(self.__sphere_list) :
            restriction = restriction & sph.restrict(initial_surface)
        return initial_surface[restriction]


class density_surface :
    DISTANCE_CUTOFF = 2.0
    RESOLUTION = 0.01
    def __init__(self, wfn, phi_res, theta_res, additional_coordinates = None, restrict = True):
        """
        :param wfn: a wavefunction for the system
        :type wfn: Wavefunction
        :param phi_res: number of points sampled at maximum contour along z axis
        :type phi_res: int
        :param theta_res: number of parallels traced
        :type theta_res: int
        :param additional_coordinates: add intermediate points to enhance sampling in certain areas.
        :type additional_coordinates: np.ndarray
        :param restrict: whether apply voronoi restriction or not
        :type restrict: bool

        """
        nuclei_coordinates, labels = wfn.getCoordinates()
        self.__wfn=wfn
        self.__phi_res=phi_res
        self.__theta_res=theta_res
        if additional_coordinates is None:
            self.__atoms_coordinates=nuclei_coordinates
        else:
            self.__atoms_coordinates = np.concatenate((nuclei_coordinates, additional_coordinates), axis=0)
        self.__apply_restrict = restrict

    def __get_number_radial_vectors(self):
        d_theta = np.pi / self.__theta_res
        phi_bin = np.array([int(np.sin(th)*self.__phi_res) for th in np.arange(self.__theta_res + 1)*d_theta])
        phi_bin[0] = 1
        phi_bin[-1] = 1
        return phi_bin


    def __create_radial_vectors(self):
        phi_bin = self.__get_number_radial_vectors()
        v = np.zeros((phi_bin.sum()-1, 3), dtype='float64')
        d_theta = np.pi / self.__theta_res
        d_phi   = 2*np.pi / phi_bin
        k = 0
        for i, th in enumerate(np.arange(self.__theta_res, dtype='float64')*d_theta):
            if phi_bin[i] == 0:
                v[k,:] = 0, 0, np.cos(th)
                k += 1
            for j, ph in enumerate(np.arange(phi_bin[i], dtype='float64')*d_phi[i]) :
                v[k, :] = np.cos(ph)*np.sin(th), np.sin(ph)*np.sin(th), np.cos(th)
                k += 1
        return v

    def restrict(self, coords, atom):
        D_ex = np.zeros((coords.shape[0], self.__atoms_coordinates.shape[0]-1), dtype='float64')
        D_in = np.linalg.norm(coords-self.__atoms_coordinates[atom, :], axis=1)
        for k, j in enumerate([i for i in range(self.__atoms_coordinates.shape[0]) if i != atom]):
            D_ex[:,k] = np.linalg.norm(coords-self.__atoms_coordinates[j,:], axis=1)

        return np.array([np.any(D_ex[i,:] < D_in[i]) for i in range(coords.shape[0])])
    def calculate_radius(self, atom, density_cont, log_tolerance = 1e-1):
        x0 = self.__atoms_coordinates[atom, :]
        v  = self.__create_radial_vectors()
        s  = np.zeros((v.shape[0],1), dtype='float64')
        h  = 1e-3
        r = np.ones((v.shape[0], 1), dtype='float64') * h
        x  = x0 + r * v
        x_up = x + h*v
        x_dw = x - h*v
        log_density_cont = np.log(density_cont)
        for i in range(10):
            p_up = np.log(self.__wfn.calculateDensity(x_up)).reshape(v.shape[0],1)
            p_dw = np.log(self.__wfn.calculateDensity(x_dw)).reshape(v.shape[0],1)
            p    = np.log(self.__wfn.calculateDensity(x)).reshape(v.shape[0],1)
            s    = np.abs(p - log_density_cont) < log_tolerance
            if np.all(s) :
                break
            der  = (p_up - p_dw)/(2*h)
            dr = (p-log_density_cont)/der
            dr[dr > 20.0] = 20.0
            dr[dr < -20.0] = -20.0
            r   -= dr
            x    = x0 + r*v
            x_up = x + h*v
            x_dw = x - h*v
        return x0 + r * v


if __name__ == "__main__" :
    # import matplotlib.pyplot as plt
    # print("testing spheres")
    # foo = np.array([
    #     [0.00000000, 0.18083455, 0.00000000],
    #     [1.56668903, -0.72251570, 0.00000000],
    #     [-1.56668903, -0.72416068, 0.00000000]
    # ], dtype='float64')
    # radius_sample = [1.0, 1.00, 1.00]
    # water_surf = surface(center_coordinates=foo, radius_list=radius_sample)
    #
    # xyz = water_surf.calculate_surface_points(phi_resolution=50, theta_resolution=50)
    #
    # sliced = np.abs(xyz[:,2]) < 0.1
    # fig, ax = plt.subplots(1)
    # ax.scatter(xyz[sliced,0], xyz[sliced,1])
    # ax.set_aspect('equal')
    # ax.set_xlabel('$\AA$')
    # ax.set_ylabel('$\AA$')
    # plt.tight_layout()
    # plt.show()
    import matplotlib.pyplot as plt

    water_wfn = Wavefunction(file=r'C:/Users/Bruno/Dropbox/doctorado/sep/input/gaussian_opt/water_gas.wfn')
    water_wfn.read()
    water_wfn.setDensityMatrix()

    water_coords, water_labels = water_wfn.getCoordinates()
    new_coords = np.zeros((2,3), dtype='float64')
    new_coords[0, :] = (water_coords[0, :] + water_coords[1,: ]) / 2.0
    new_coords[1, :] = (water_coords[0, :] + water_coords[2, :]) / 2.0
    ddd = density_surface(water_wfn, phi_res=40, theta_res=20, additional_coordinates=new_coords)
    xyz = ddd.calculate_radius(atom=0, density_cont=0.1)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], color = 'blue')
    xyz = ddd.calculate_radius(atom=1, density_cont=0.1)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], color='red')
    xyz = ddd.calculate_radius(atom=2, density_cont=0.1)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], color='red')
    xyz = ddd.calculate_radius(atom=3, density_cont=0.1)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], color='purple')
    xyz = ddd.calculate_radius(atom=4, density_cont=0.1)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], color='purple')
    plt.show()

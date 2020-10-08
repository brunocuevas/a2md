import time
import sys
import numpy as np
from scipy import signal as sig


def create_orthogonal_basis(main_axis):
    main_axis /= np.linalg.norm(main_axis)
    second_axis = np.array(
        [
            -main_axis[1],
            main_axis[0],
            0
        ]
    )
    third_axis = np.cross(main_axis, second_axis)
    second_axis /= np.linalg.norm(second_axis)
    third_axis /= np.linalg.norm(third_axis)
    basis = np.array([main_axis, second_axis, third_axis])
    return basis.T


def get_angle(v1, v2):
    m1 = np.linalg.norm(v1)
    m2 = np.linalg.norm(v2)
    tt = np.arccos((np.sum(v1 * v2)) / (m1 * m2))
    return tt


def get_projected_angle(v1, v2, axis):
    u1 = np.array([v1[i] for i in range(3) if i != axis])
    u2 = np.array([v2[i] for i in range(3) if i != axis])
    v1_0 = v1.copy()
    v2_0 = v2.copy()
    v1_0[axis] = 0
    v2_0[axis] = 0
    sign = np.cross(v1_0, v2_0)[axis]
    if sign > 0:
        return get_angle(u1, u2)
    else:
        return -get_angle(u1, u2)


def rotz(coords, angle):
    rz = np.array(
        [[np.cos(angle), -np.sin(angle), 0],
         [np.sin(angle), np.cos(angle), 0],
         [0, 0, 1]],
    )
    return rz.dot(coords)


def rotx(coords, angle):
    rx = np.array(
        [[1, 0, 0],
         [0, np.cos(angle), -np.sin(angle)],
         [0, np.sin(angle), np.cos(angle)]]
    )
    return rx.dot(coords)


def roty(coords, angle):
    ry = np.array(
        [[np.cos(angle), 0, np.sin(angle)],
         [0, 1, 0],
         [-np.sin(angle), 0, np.cos(angle)]]
    )
    return ry.dot(coords)


def write_pdb(filename, coordinates):
    format_pdb = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
    with open(filename, 'w') as f:
        for i in range(coordinates.shape[0]):
            f.write(
                format_pdb.format(
                    'HETATM', i, 'He', '',
                    '', '', 0, '',
                    coordinates[i, 0], coordinates[i, 1], coordinates[i, 2],
                    0, 0, 'He', ''
                )
            )


# pyProb defined Errors

class Error(Exception):
    pass


class ConvolutionError(Error):
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


# Classes

class VolumeBaseClass:
    def __init__(self, name, verbose=True, file=None):
        self.__name = name
        self.__verbose = verbose
        self.__file = file
        self.__f = None
        self.__flags = dict(
            log="@",
            dat="#"
        )

    def log(self, message):
        ltime = time.localtime()
        if self.__verbose:
            print(
                "@ [%s] %s - %02d:%02d:%02d" %
                (self.__name, message, ltime.tm_hour, ltime.tm_min, ltime.tm_sec)
            )

    def dat(self, message):
        print("# [%s] %s" % (self.__name, message))

    def warning(self, message):
        ltime = time.localtime()
        print(
            "$ [%s] %s - %02d:%02d:%02d" % (
                self.__name, message, ltime.tm_hour, ltime.tm_min, ltime.tm_sec
            )
        )

    def open_log_file(self):
        self.__f = open(self.__file, 'a')

    def log2file(self, messsage, kind="log"):
        if self.__f is None:
            pass
        else:
            string2write = "%s[%s] %s\n" % (self.__flags[kind], self.__name, messsage)
            self.__f.write(string2write)


class Volume(VolumeBaseClass):
    def __init__(
            self,
            filename=None,  # filename of dx file that can be read
            dxvalues=None,  # numpy 3D tensor
            r0=None,  # system origin of coordinates
            verbose=True,  # level of detail in the logfile
            basis=None,  # system coordinates basis
            masscenter=None,
            phi=None,
            psi=None,
            edge=0
    ):
        VolumeBaseClass.__init__(self, 'volumetric', verbose=verbose)
        # The object can be defined by reading a tensor defined in python
        # or by reading a volumetric dx file. In case that it is read from a
        # tensor, it automatically sets the attributes according to the
        # features of that tensor. Otherwise, they are set to None until
        # file is read.
        self.__fn = filename
        self.__dx = dxvalues
        if dxvalues is None:
            self.__nx = None
            self.__ny = None
            self.__nz = None
        else:
            self.__nx, self.__ny, self.__nz = dxvalues.shape
        if r0 is None and not (dxvalues is None):
            self.__r0 = np.array([0, 0, 0])
            self.log('setting r0 coordinates by default')
        elif r0 is None and not (filename is None):
            self.__r0 = None
        else:
            self.__r0 = r0
        if basis is None and not (dxvalues is None):
            self.__X = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            self.log('setting coordinates basis by default')
        elif basis is None and not (filename is None):
            self.__X = None
        else:
            self.__X = basis
        if dxvalues is None:
            self.shape = (0, 0, 0)
        else:
            self.shape = dxvalues.shape
        self.__lb = None
        self.__ub = None
        self.__v = verbose
        if phi is None:
            self.__phi = None
        else:
            self.__phi = phi
        if psi is None:
            self.__psi = None
        else:
            self.__psi = psi
        if masscenter is None:
            self.__massCenter = None
        else:
            self.__massCenter = masscenter
        self.__edge = edge

    # Default Methods
    def __iter__(self):
        return self

    def __next__(self):
        try:
            lb = self.__lb[0]
            ub = self.__ub[0]
            self.log("lb : %d ub : %d" % (lb, ub))
        except IndexError:
            raise StopIteration
        try:
            self.__lb = self.__lb[1:]
            self.__ub = self.__ub[1:]
        except IndexError:
            self.__lb = []
            self.__ub = []

        x_size = int(ub - lb)
        n = self.shape[0]
        y_size = self.shape[1]
        z_size = self.shape[2]
        dx = np.zeros((x_size, y_size, z_size))
        edge = 0
        if ub < 0:
            edge -= 1
        elif lb < 0 < ub:
            edge -= 1
            dx[-ub:, :, :] = self.__dx[:ub, :, :]
        elif lb > 0 and ub < n:
            dx[:, :, :] = self.__dx[lb:ub, :, :]
        elif lb < n < ub:
            edge += 1
            nlb = n - lb
            dx[:nlb, :, :] = self.__dx[lb:, :, :]
        elif lb > n:
            pass
        origin = self.__r0 + self.__X.dot(lb * np.array([1, 0, 0]))
        return Volume(
            dxvalues=dx,
            r0=origin,
            basis=self.__X.copy(),
            verbose=self.__v,
            masscenter=self.__massCenter,
            phi=self.__phi,
            psi=self.__psi,
            edge=edge
        )

    def __str__(self):
        if (self.__fn is None) and (self.__dx is None):
            return "volumetric "
        elif self.__dx is None:
            return "volumetric file : %s" % self.__fn
        elif self.__fn is None:
            return "volumetric shape : %d %d %d" % (
                self.shape[0], self.shape[1], self.shape[2]
            )
        else:
            return "volumetric file : %s shape : %d %d %d" % (
                self.__fn, self.shape[0], self.shape[1], self.shape[2]
            )

    # Static Methods
    @staticmethod
    def cos(angle):
        return np.cos(-angle)

    @staticmethod
    def sin(angle):
        return np.sin(-angle)

    # Private Methods
    def __get_center(self):

        total = np.sum(self.__dx)
        cx = 0
        cy = 0
        cz = 0
        for ix in range(self.__nx):
            for iy in range(self.__ny):
                for iz in range(self.__nz):
                    ratio = self.__dx[ix, iy, iz] / total
                    cx += ratio * ix
                    cy += ratio * iy
                    cz += ratio * iz

        self.__massCenter = np.array([cx, cy, cz]) + self.__r0
        return np.array([cx, cy, cz])

    def __interpolate(self, coordinates, interpolate=True):
        ndx = np.zeros((self.__nx, self.__ny, self.__nz))
        ix = 0;
        iy = 0;
        iz = 0
        n = 0
        r0 = self.get_r0()
        for c in coordinates:
            cx = int(c[0] - r0[0])
            cy = int(c[1] - r0[1])
            cz = int(c[2] - r0[2])
            if any([(cx >= (self.__nx - 1)), (cy >= (self.__ny - 1)), (cz >= (self.__nz - 1))]):
                ndx[ix, iy, iz] = 0.0
                n += 1
            elif any([(cx < 1), (cy < 1), (cz < 1)]):
                ndx[ix, iy, iz] = 0.0
                n += 1
            else:
                deltax = c[0] - float(cx)
                deltay = c[1] - float(cy)
                deltaz = c[2] - float(cz)
                ndx[ix, iy, iz] = self.__dx[cx, cy, cz]
                if interpolate:
                    ndx[ix, iy, iz] = (ndx[ix, iy, iz] + deltax * self.__dx[cx + 1, cy, cz]) / (1 + deltax)
                    ndx[ix, iy, iz] = (ndx[ix, iy, iz] + deltay * self.__dx[cx, cy + 1, cz]) / (1 + deltay)
                    ndx[ix, iy, iz] = (ndx[ix, iy, iz] + deltaz * self.__dx[cx, cy, cz + 1]) / (1 + deltaz)
            iz += 1
            if iz == self.__nz:
                iy += 1
                iz = 0
            if iy == self.__ny:
                iy = 0
                ix += 1
        return ndx

    def __rotx(self, coordinates, angle):
        rx = np.array(
            [[1, 0, 0],
             [0, self.cos(angle), -self.sin(angle)],
             [0, self.sin(angle), self.cos(angle)]]
        )
        return coordinates.dot(rx.T)

    def __roty(self, coordinates, angle):
        ry = np.array(
            [[self.cos(angle), 0, self.sin(angle)],
             [0, 1, 0],
             [-self.sin(angle), 0, self.cos(angle)]]
        )
        return coordinates.dot(ry.T)

    def __rotz(self, coordinates, angle):
        rz = np.array(
            [[self.cos(angle), -self.sin(angle), 0],
             [self.sin(angle), self.cos(angle), 0],
             [0, 0, 1]],
        )
        return coordinates.dot(rz.T)

    def __set_coords_4_rotation(self):
        i = 0
        coordinates = np.zeros((self.__nx * self.__ny * self.__nz, 3))
        r0 = self.get_r0()
        for ix in range(self.__nx):
            for iy in range(self.__ny):
                for iz in range(self.__nz):
                    coordinates[i, 0] = float(ix)
                    coordinates[i, 1] = float(iy)
                    coordinates[i, 2] = float(iz)
                    coordinates[i, :] += r0
                    i += 1
        return coordinates

    # Public Methods
    def add(self, volume):
        self.log('adding tensor')
        try:
            cnter = volume.get_geometric_center()
            shape = volume.shape
            volume = volume.get_volume()
        except AttributeError:
            raise IOError("not a volumetric instance")

        cx = int((cnter[0] - self.__r0[0]) / self.__X[0, 0])
        cy = int((cnter[1] - self.__r0[1]) / self.__X[1, 1])
        cz = int((cnter[2] - self.__r0[2]) / self.__X[2, 2])
        r = np.array(shape) % 2

        d = np.array([
            cx - int(shape[0] / 2),
            cy - int(shape[1] / 2),
            cz - int(shape[2] / 2)
        ], dtype='int64')

        e = np.array([
            cx + int(shape[0] / 2),
            cy + int(shape[1] / 2),
            cz + int(shape[2] / 2)
        ], dtype='int64') + r

        d[d < 0] = 0
        e[e < 0] = 0
        n = np.array([self.__nx, self.__ny, self.__nz])
        d[d > n] = n[d > n]
        e[e > n] = n[e > n]

        try:
            self.__dx[
            d[0]: e[0],
            d[1]: e[1],
            d[2]: e[2]
            ] += volume
        except IndexError:
            raise IOError("volumes should have the same size")

    def convolve(self, movingelement):
        """
        convolve
        takes a tensor or a volumetric file, and performs a discrete convolution.
        :param movingelement:
        :return: convolved tensor, coordinates of the maximum of the convolution
        """
        self.log("starting convolution in real space")

        try:
            dx = self.__dx
        except AttributeError:
            raise AttributeError("must read a volumetric file before convolution")
        try:
            movingelement_tnsor = movingelement.get_volume()
            self.log('tensor obtained from volumetric instance')
        except AttributeError:
            movingelement_tnsor = movingelement

        if movingelement_tnsor.shape[0] > self.shape[0]:
            size_dif = movingelement_tnsor.shape[0] - self.shape[0]
            movingelement_tnsor[- size_dif - 1, :, :] += np.sum(movingelement_tnsor[-size_dif:, :, :], axis=0)
            movingelement_tnsor = movingelement_tnsor[: -size_dif, :, :]
        if movingelement_tnsor.shape[0] > self.shape[0]:
            print("ERROR : mets : %d self : %d" % (movingelement_tnsor.shape[0], self.shape[0]))

        # AVOIDING ZERO-PADDING - The result of convolution should have only 2D
        cnvlvd = sig.correlate(dx, movingelement_tnsor, method='auto', mode='valid')[-1, :, :]
        cnvlvd /= (movingelement_tnsor.shape[0] * movingelement_tnsor.shape[1] * movingelement_tnsor.shape[2])
        if len(cnvlvd.shape) > 2:
            raise ConvolutionError(
                expression="cnvlvd /= (movingelement_tnsor.shape[0] * movingelement_tnsor.shape[1] * movingelement_tnsor.shape[2])",
                message="resulting tensor of the convolution has 3-dimensional shape"
            )
        rny = np.arange(0, self.shape[1] - movingelement_tnsor.shape[1] + 1, dtype='float64')
        rnz = np.arange(0, self.shape[2] - movingelement_tnsor.shape[2] + 1, dtype='float64')
        try:
            e = np.where(cnvlvd == cnvlvd.max())
        except ValueError:
            raise RuntimeError
        if len(e[0]) > 1:
            rndm = np.random.randint(0, high=len(e[0]))
            cy = e[0][rndm]
            cz = e[1][rndm]
        else:
            cy, cz = np.unravel_index(cnvlvd.argmax(), cnvlvd.shape)

        shift_y = rny[cy]
        shift_z = rnz[cz]

        r0_x, r0_y, r0_z = self.get_edge_coords()
        r = np.array(
            [
                r0_x,
                r0_y + shift_y,
                r0_z + shift_z
            ],
            dtype='float64'
        )
        topconv_x = r[0]
        topconv_y = r[1]
        topconv_z = r[2]

        self.log(
            "convolution done. the maximum was found at %f %f" %
            (
                topconv_y,
                topconv_z,
            )
        )
        return cnvlvd[cy, cz], topconv_x, topconv_y, topconv_z

    ## Getters
    def get_basis(self):
        if self.__X is None:
            return np.identity(3)
        else:
            return self.__X.copy()

    def get_coordinates_under_treshold(self, treshold=0.5):

        dx = self.__dx.copy()
        nvals = dx[dx > treshold].size
        r = np.zeros((nvals, 3), dtype='float64')
        k = 0
        for ix in range(self.shape[0]):
            for iy in range(self.shape[1]):
                for iz in range(self.shape[2]):
                    if dx[ix, iy, iz] > treshold:
                        r[k, :] = self.__X.dot(np.array([ix, iy, iz], dtype='float64')) + self.__r0
                        k += 1

        return r

    def get_difussion_axis(self, treshold=0.5):
        """

        :param treshold:
        :return:
        """
        r = self.get_coordinates_under_treshold(treshold=treshold)
        r = r - r.mean(axis=0)
        c = (1 / r.shape[0]) * r.T.dot(r)
        evl, evc = np.linalg.eig(c)

        return evl, evc

    def get_edge_coords(self):
        """
        get Edge Coordinates
        returns the coordinates of the upper corner of the tensor
        :return: x_edge, y_edge, z_edge
        """
        r0_x = self.__r0[0]  # - ((self.__nx / 2) * self.__X[0, 0])
        r0_y = self.__r0[1]  # - ((self.__ny / 2) * self.__X[1, 1])
        r0_z = self.__r0[2]  # - ((self.__nz / 2) * self.__X[2, 2])
        return r0_x, r0_y, r0_z

    def get_geometric_center(self):
        return np.array([
            self.__nx / 2,
            self.__ny / 2,
            self.__nz / 2
        ]).dot(self.__X) + self.__r0

    def get_mass_center(self):
        return self.__get_center()

    def get_r0(self):
        return self.__r0.copy()

    def get_slice_number(self):
        return len(self.__lb)

    def get_volume(self):
        """
        getVolume
        returns a volumetric file obtained by slicing the tensor
        :return: volumetric
        """
        return self.__dx

    def multiply(self, factor):
        self.__dx *= factor

    def eval(self, fun):
        res = self.__X[0, 0]
        minx, miny, minz = self.__r0
        xx = minx + (np.arange(self.__nx) * res)
        yy = miny + (np.arange(self.__ny) * res)
        zz = minz + (np.arange(self.__nz) * res)

        # dx = np.zeros((xx.size, yy.size, zz.size))

        X, Y, Z = np.meshgrid(xx, yy, zz)
        r = np.stack([Z.flatten(), X.flatten(), Y.flatten()], axis=1)
        dx = fun(r).reshape(xx.size, yy.size, zz.size).T
        self.__dx = dx

    def read(self):
        """
        read
        reads a dx file and sets the tensor and the other attributes
        :return:
        """
        if self.__fn is None:
            raise IOError("filename is not defined")
        self.log("reading %s" % self.__fn)
        dx = None
        ix = 0
        iy = 0
        iz = 0
        i = 0
        x = np.zeros((3, 3))
        ny = None
        nz = None
        try:
            f = open(self.__fn)
        except FileNotFoundError:
            self.log("file %s was not found. dying" % self.__fn)
            sys.exit(1)
        for line in f:
            if line[0] == '#':
                continue
            terms = line.split()
            if terms[0] == 'object':
                if terms[3] == 'gridpositions':
                    self.__nx = int(terms[5])
                    self.__ny = int(terms[6])
                    self.__nz = int(terms[7])
                    nx = int(terms[5])
                    ny = int(terms[6])
                    nz = int(terms[7])
                    dx = np.zeros((
                        int(terms[5]),
                        int(terms[6]),
                        int(terms[7]),
                    ))
                    self.log(
                        'volumetric tensor was created. Size x %d y %d z %d' % (
                            nx, ny, nz
                        )
                    )
                else:
                    pass
            elif terms[0] == 'origin':
                self.__r0 = np.array(
                    [float(terms[1]), float(terms[2]), float(terms[3])]
                )
                self.log(
                    'origin coordinates were set to %4.3f %4.3f %4.3f' % (
                        float(terms[1]), float(terms[2]), float(terms[3])
                    )
                )
            elif terms[0] == 'delta':
                x[i, i] = terms[i + 1]
                i += 1
            else:
                for t in terms:
                    try:
                        dx[ix, iy, iz] = float(t)
                    except TypeError:
                        raise IOError('dx file format was not correct')
                    except ValueError:
                        break
                    iz += 1
                    if iz == nz:
                        iy += 1
                        iz = 0
                    if iy == ny:
                        ix += 1
                        iy = 0
        self.__dx = dx
        self.shape = dx.shape
        self.__X = x
        self.log("dx reading is finished")
        return True

    def rotate_x(self, coordinates, angle):
        return self.__rotx(coordinates, angle)

    def rotate_y(self, coordinates, angle):
        return self.__roty(coordinates, angle)

    def rotate_z(self, coordinates, angle):
        return self.__rotz(coordinates, angle)

    def rotate(self, newbasis, interpolate=True):
        """
        rotate allows the rotation around the --mass center-- of the volume,
        to align an arbitrary axis x with the x axis of the tensor. If interpolate
        is set to True, it performs a linear approximation of the values at each given
        point
        :param newbasis:
        :param interpolate:
        :return:
        """
        self.log('changing basis system')

        ## Rotation of tensor coordinates
        coordinates = self.__set_coords_4_rotation()
        coordinates = coordinates.dot(newbasis.T)
        self.log('interpolating volume tensor')
        ndx = self.__interpolate(coordinates, interpolate=interpolate)
        self.__dx = ndx
        self.log('interpolation, done')
        return True

    def rotate_around_x(self, angle, interpolate=True):
        coordinates = self.__set_coords_4_rotation()
        coordinates = self.__rotx(coordinates, angle)
        ndx = self.__interpolate(coordinates, interpolate=interpolate)
        self.__dx = ndx

    def set_basis(self, basis):
        self.__X = basis

    def set_mass_center(self, mass_center):
        self.__massCenter = mass_center
        return self.__massCenter

    def set_r0(self, r0):
        self.__r0 = r0

    def set_volume(self, dx):
        self.__dx = dx
        self.__nx, self.__ny, self.__nz = dx.shape
        self.shape = dx.shape
        if self.__X is None:
            self.__X = np.identity(3, dtype='float64')
        if self.__r0 is None:
            self.__r0 = np.zeros(3)

    def slice(self, application, step, extend, overlap):
        """
        slicing
        allows to cut the tensor in overlapping boxes, so they
        can be obtained by iterating over the object
        :param application: point where the slicing starts
        :param step: distance of the starting point of slicing
        :param extend: total distance sliced
        :param overlap: extension of each slice
        :return:
        """
        self.log("starting slicing")
        deltax = step
        n_steps = int((extend - application) / deltax) + 1
        # start   = application
        start = np.round(application / np.linalg.norm(self.__X[0, :]), decimals=0).astype('int')
        low_bounds = [start]
        up_bounds = [start + overlap]
        overlap = overlap * self.__X[0, 0]
        for i in range(n_steps):
            start = start + step
            end = start + overlap
            low_bounds.append(int(start))
            up_bounds.append(int(end))
        self.__lb = low_bounds
        self.__ub = up_bounds
        return True

    def write(self, filename):
        f = open(filename, 'w')
        f.write('#\n')
        f.write(
            'object 1 class gridpositions counts %d %d %d\n' % (self.__nx, self.__ny, self.__nz)
        )
        f.write('origin %f %f %f\n' % (self.__r0[0], self.__r0[1], self.__r0[2]))
        f.write('delta %5.4f %5.4f %5.4f\n' % (self.__X[0, 0], self.__X[0, 1], self.__X[0, 2]))
        f.write('delta %5.4f %5.4f %5.4f\n' % (self.__X[1, 0], self.__X[1, 1], self.__X[1, 2]))
        f.write('delta %5.4f %5.4f %5.4f\n' % (self.__X[2, 0], self.__X[2, 1], self.__X[2, 2]))
        f.write('object 2 class gridconnections counts %d %d %d\n' % (
            self.__nx, self.__ny, self.__nz
        )
                )
        f.write('object 3 class array type double rank 0 items %d data follows\n' % (self.__nx * self.__ny * self.__nz))
        ll = []
        for ix in range(self.__nx):
            for iy in range(self.__ny):
                for iz in range(self.__nz):
                    ll.append(self.__dx[ix, iy, iz])
                    if len(ll) == 3:
                        f.write('{:12.11e} {:12.11e} {:12.11e}'.format(ll[0], ll[1], ll[2]) + '\n')
                        ll = []
        if len(ll) > 0:
            try:
                f.write('{:12.11e} {:12.11e}'.format(ll[0], ll[1]) + '\n')
            except IndexError:
                f.write('{:12.11e}'.format(ll[0]) + '\n')
        f.close()

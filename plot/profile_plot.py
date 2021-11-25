## generate contourline from 2d-data
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
import math
from functools import singledispatch

@singledispatch
def linear_interpolate_1d(xl_, xr_, t_):
    return xl_ * (1.0 - t_) + xr_ * t_


@linear_interpolate_1d.register(float)
@linear_interpolate_1d.register(np.float64)
@linear_interpolate_1d.register(np.float32)
def _(xl_: float, xr_: float, t_: float) -> float:
    return xl_ * (1.0 - t_) + xr_ * t_


@linear_interpolate_1d.register(np.ndarray)
def _(xl_: np.ndarray, t_: float) -> float:
    xi = np.min([math.floor(t_),xl_.shape[0]-1])
    xj = xi + 1
    xt = t_ - xi
    return linear_interpolate_1d(xl_[xi], xl_[xj], xt)


def linear_interpolate_2d(xl_: np.ndarray, xt_: float, yt_: float) -> float:
    xi = np.min([math.floor(xt_), xl_.shape[0]-2])
    yi = np.min([math.floor(yt_), xl_.shape[1]-2])
    x_tmp_0 = np.array([[xl_[i + xi, j + yi] for i in [0, 1]] for j in [0, 1]])
    x_tmp_1 = [linear_interpolate_1d(x_tmp_0[0, i], x_tmp_0[1, i], xt_ - xi) for i in [0, 1]]
    return linear_interpolate_1d(x_tmp_1[0], x_tmp_1[1], yt_ - yi)


def linear_interpolate_2ds(xl_: np.ndarray, xt_: np.ndarray, yt_: np.ndarray) -> np.ndarray:
    x_rt = np.ones_like(xt_)
    for k,tmp in enumerate(x_rt):
        x_rt[k] = linear_interpolate_2d(xl_,xt_[k],yt_[k])
    return x_rt


class read_data:
    def __init__(self, path_, step_, tag_=2, rx_=None, ry_=None, rz_=None):
        if rz_ is None:
            rz_ = [None, None]
        if ry_ is None:
            ry_ = [None, None]
        if rx_ is None:
            rx_ = [None, None]
        self.step_ = step_
        self.data_path_ = "%s/Plane.data.step.%d." % (path_, step_)
        self.coor_path_ = "%s/Plane.coord.step.%d." % (path_, step_)
        self.tag_ = tag_
        self._sxl, self._sxr = rx_
        self._syl, self._syr = ry_
        self._szl, self._szr = rz_

        self.Xx = None
        self.Xy = None
        self.Xz = None
        self.Xd = None
        self.XC = None
        self.XD = None

        self.__other = {}

    def load(self):

        self.XC = np.loadtxt(self.coor_path_ + str(self.tag_), delimiter=",", dtype=np.string_, skiprows=1)
        self.XD = np.loadtxt(self.data_path_ + str(self.tag_), delimiter=",", dtype=np.string_)
        X = (self.XD[:, :-1]).astype(np.float32)
        x_mask = self.XC[:, 3].reshape(X.shape).astype(np.int32)
        self.Xx = self.XC[:, 0].reshape(X.shape).astype(np.float32)
        self.Xy = self.XC[:, 1].reshape(X.shape).astype(np.float32)
        self.Xz = self.XC[:, 2].reshape(X.shape).astype(np.float32)
        self.Xd = np.ma.array(X, mask=(x_mask < 0))
        plot_x = None
        plot_y = None
        plot_d = None

        if (self.tag_ in [0, 4]) and (self._sxl is not None):
            self._sxl, self._sxr = np.searchsorted(self.Xx[0, :], [self._sxl, self._sxr])
            self._szl, self._szr = np.searchsorted(self.Xz[:, 0], [self._szl, self._szr])
            plot_x = self.Xx[self._szl:self._szr, self._sxl:self._sxr]
            plot_y = self.Xz[self._szl:self._szr, self._sxl:self._sxr]
            plot_d = self.Xd[self._szl:self._szr, self._sxl:self._sxr]
        elif (self.tag_ in [0, 4]) and (self._sxl is None):
            plot_x = self.Xx
            plot_y = self.Xz
            plot_d = self.Xd
        elif (self.tag_ in [1, 2, 3]) and (self._sxl is not None):
            self._sxl, self._sxr = np.searchsorted(self.Xx[:, 0], [self._sxl, self._sxr])
            self._syl, self._syr = np.searchsorted(self.Xy[0, :], [self._syl, self._syr])
            plot_x = self.Xx[self._sxl:self._sxr, self._syl:self._syr]
            plot_y = self.Xy[self._sxl:self._sxr, self._syl:self._syr]
            plot_d = self.Xd[self._sxl:self._sxr, self._syl:self._syr]
        elif (self.tag_ in [1, 2, 3]) and (self._sxl is None):
            plot_x = self.Xx
            plot_y = self.Xy
            plot_d = self.Xd

        self.__other['plot'] = [plot_x, plot_y, plot_d]

    def __del__(self):
        del self.XC
        del self.XD

    # noinspection PyShadowingNames
    def plot(self, ax_):

        ax_.set_xlabel("X[m]")
        ax_.set_ylabel("Y[m]")

        plot_x, plot_y, plot_d = self.__other['plot']

        cx = ax_.contourf(plot_x, plot_y, plot_d, cmap="bwr", levels=30)
        px = plt.colorbar(cx, fraction=0.025)

        if self.tag_ == 2:
            # height profile
            # ax_.contour(plot_x, plot_y, plot_d, levels=[-800.0, 0.0, 800.0])
            px.set_label("Height[m]")
        elif self.tag_ == 0:
            # XZ-density
            # print(self.Xx[0, :], np.diff(self.Xx[0, :]))
            px.set_label("Density")
        elif self.tag_ == 1:
            # XY-density
            px.set_label("Density")

        return px

    def plot_3d(self, ax_):

        ax_.set_xlabel("X[m]")
        ax_.set_ylabel("Y[m]")

        plot_x, plot_y, plot_d = self.__other['plot']

        cx = ax_.plot_surface(plot_x, plot_y, plot_d, cmap="bwr")
        px = plt.colorbar(cx)

        if self.tag_ == 2:
            # height profile
            # ax_.contour(plot_x, plot_y, plot_d, levels=[-800.0, 0.0, 800.0])
            px.set_label("Height[m]")
        elif self.tag_ == 0:
            # XZ-density
            # print(self.Xx[0, :], np.diff(self.Xx[0, :]))
            px.set_label("Density")
        elif self.tag_ == 1:
            # XY-density
            px.set_label("Density")

        return px

    def get_plot(self):
        return self.__other['plot']


class cngn:
    data = None
    x, y, z = [None, None, None]
    segments = None
    contours = []
    level = None

    def __init__(self, x_=None, y_=None, z_=None):
        self.x, self.y, self.z = x_, y_, z_
        self.segments = None if z_ is None else measure.find_contours(z_, 0)

    def reinit(self, other: read_data, level=0.):
        self.x, self.y, self.z = other.get_plot()
        self.level = level
        self.segments = None if self.z is None else measure.find_contours(self.z, level)
        # print(type(self.segments))
        return self.segments

    def contour(self, level: float):
        self.level = level
        self.segments = None if self.z is None else measure.find_contours(self.z, level)
        return self.segments

    def show(self, ax: plt.axis) -> plt.axis:
        ax.imshow(self.z.T, origin='lower')
        for cns in self.segments:
            ax.plot(cns[:, 0], cns[:, 1])
        return ax

    def plot_lines(self, ax: plt.axis, indices: list = None, colors="k") -> plt.axis:
        if indices is None:
            for cns in self.contours:
                ax.plot(cns[0], cns[1])
        else:
            for index in indices:
                ax.plot(self.contours[index][0], self.contours[index][1], colors)

        return ax

    def contour_lines(self):
        self.contours = []
        for cns in self.segments:
            self.contours.append([linear_interpolate_2ds(self.x, cns[:, 0], cns[:, 1]),
                                  linear_interpolate_2ds(self.y, cns[:, 0], cns[:, 1])])

    def __call__(self, x_: float, y_: float) -> list:
        return [linear_interpolate_2d(self.x, x_, y_), linear_interpolate_2d(self.y, x_, y_)]

    def find_nearest(self, px_:float, py_:float):
        distances = []
        for cns in self.contours:
            tx1_, tx2_ = cns
            td0_ = np.sqrt((tx1_ - px_)**2 + (tx2_ - py_)**2)
            distances.append([np.argmax(td0_), np.max(td0_), np.argmin(td0_), np.min(td0_), np.mean(td0_)])

        distances = np.array(distances)
        # print(distances)
        # exit(0)
        return np.argmin(distances[:, 3])

    def get_extent(self, index: int) -> list:
        if -1 == index:
            return [self.z.min(), self.z.max()]
        tmp_x = self.contours[index][0]
        tmp_y = self.contours[index][1]
        Lx = np.max(tmp_x) - np.min(tmp_x)
        Ly = np.max(tmp_y) - np.min(tmp_y)
        Depth = self.z.min()
        return [Lx, Ly, self.level - Depth, Lx / Ly, Ly / Lx]



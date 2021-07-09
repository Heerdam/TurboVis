"""The WaveBlocks Project

Function for mapnp.ping complex numbers to colors specified
by the usual color map used in quantum mechanics.

@author: R. Bourquin
@copyright: Copyright (C) 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb


def color_map(data, phase=None, modulus=None, darken=1.0):
    """Color mapnp.ping according to the QM standard map.

    :param data: The complex numbers.
    :param phase: The phase of the complex numbers, computed if not given.
    :param modulus: The modulus of the complex numbers, computed if not given.
    :param darken: How strong to take into account the modulus of the data to darken colors.
                   Values with :math:`|z| = R` will get fully saturated colors
                   while :math:`|z| = 0` is black and :math:`|z| \rightarrow \infty`
                   get whiter and whiter.
    :type darken: Float or ``None`` to disable darkening of colors. Default is :math:`R = 1.0`.
    :param axes: The axes instance used for plotting.
    """
    if len(data.shape) == 1:
        hsv_colors = np.empty((1,) + data.shape + (3,))
    else:
        hsv_colors = np.empty(data.shape + (3,))

    if phase is None:
        phase = np.angle(data)

    hsv_colors[..., 0] = 0.5 * np.fmod(phase + 2.0 * np.pi, 2.0 * np.pi) / np.pi
    hsv_colors[..., 1] = 1.0
    hsv_colors[..., 2] = 1.0

    # Darken colors such that 0+0i maps to black
    if darken is not None:
        if modulus is None:
            modulus = np.abs(data)

        # Lightness
        hsv_colors[..., 2] = 2.0 * np.arctan2(np.real(modulus), darken) / np.pi

        # Saturation
        l = hsv_colors[..., 2]
        hsv_colors[..., 1] = np.where(l <= 0.5, 2 * l, 2 * (1 - l))

    return hsv_to_rgb(hsv_colors)

"""The WaveBlocks Project

Function for plotting complex valued functions
of two real variables with the values encoded
by the usual color code.

@author: R. Bourquin
@copyright: Copyright (C) 2012, 2014, 2016 R. Bourquin
@license: Modified BSD License
"""
def plotcf2d(x, y, z, darken=None, axes=None, limits=None, **kwargs):
    r"""Plot complex valued functions :math:`\mathbb{R}^2 \rightarrow \mathbb{C}`
    with the usual color code.

    :param x: The :math:`x` values.
    :param x: The :math:`y` values.
    :param z: The values :math:`z = f(x,y)`.
    :param darken: How strong to take into account the modulus of the data to darken colors.
                   Values with :math:`|z| = R` will get fully saturated colors
                   while :math:`|z| = 0` is black and :math:`|z| \rightarrow \infty`
                   get whiter and whiter.
    :type darken: Float or ``None`` to disable darkening of colors. Default is :math:`R = 1.0`.
    :param axes: The axes instance used for plotting.
    """
    if limits is None:
        xmin = np.real(x).min()
        xmax = np.real(x).max()
        ymin = np.real(y).min()
        ymax = np.real(y).max()
        extent = [xmin, xmax, ymin, ymax]
    else:
        xmin = limits[0]
        xmax = limits[1]
        ymin = limits[2]
        ymax = limits[3]
        extent = [xmin, xmax, ymin, ymax]

    kw = {'extent': extent,
          'origin': 'lower',
          'interpolation': 'bilinear',
          'aspect': 'equal',
          'cmap': 'hsv',
          'vmin': 0,
          'vmax': 2.0 * np.pi}
    kw.update(kwargs)

    # Plot to the given axis instance or retrieve the current one
    if axes is None:
        axes = plt.gca()

    # Region to cut out
    x = x.reshape(1, -1)
    y = y.reshape(-1, 1)
    i = np.where((xmin <= x) & (x <= xmax))[1]
    j = np.where((ymin <= y) & (y <= ymax))[0]
    I, J = np.meshgrid(i, j)

    # Color code and plot the data
    return axes.imshow(color_map(z[J, I], darken=darken), **kw)#, cmap

def phaseplot(x, y, Z, filename, darken=0.05, bbox_inches=None):
	fig = plt.figure()
	fig.add_subplot(111)
	plt.title("Wavefunction")
	cs = plotcf2d(x, y, Z, darken=darken)
	cbar = plt.colorbar(cs, ticks=[0.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0, 2.0 * np.pi])
	cbar.ax.set_yticklabels([r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
	fig.savefig(filename, bbox_inches=bbox_inches)
	plt.close(fig)

def Wavefunction(X, Y):
    q = np.array([2.0, 1.0])    # position center
    p = np.array([-0.5, -1.0])  # momentum center
    sigma = 1.0                 # standard deviation

    d = 2                       # spacial dimension
    normalizing_factor = np.pi**(-d / 4.0) / np.sqrt(sigma)

    Xq = X - q[0]
    Yq = Y - q[1]
    phase = p[0] * Xq + p[1] * Yq

    return normalizing_factor * np.exp(-0.5 * (Xq**2 + Yq**2) / sigma + 1.0j * phase)

x = np.linspace(-6.0, 6.0, 2**12, endpoint=False)
y = np.linspace(-6.0, 6.0, 2**12, endpoint=False)
X, Y = np.array(np.meshgrid(x, y))
Z = Wavefunction(X, Y)
phaseplot(x, y, Z, 'wavefunction.png')

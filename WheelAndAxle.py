r"""
Definition
----------

Calculates WheelAndAxle.



References
----------

Authorship and Verification
---------------------------

* **Author:** --- **Date:** 2018YYY-04m-18d
* **Last Modified by:** --- **Date:** 2018YYY-04m-18d
* **Last Reviewed by:** --- **Date:** 2018YYY-04m-18d
"""

from numpy import inf
import numpy as np
from scipy.special import spherical_jn
from scipy.integrate import quad, dblquad

name = "WheelAndAxle"
title = "User model for WheelAndAxle"
description = """"""

parameters = [
#   ["name", "units", default, [lower, upper], "type", "description"],
    ['AxleRadius', 'Angstrom', 80.0, [0, inf], '', ''],
    ['AxleLength', 'Angstrom', 1.0, [0, inf], '', ''],
    ['AxleSld', '', 1.0, [-inf, inf], '', ''],
    ['WheelRadius', 'Angstrom', 1.0, [0, inf], '', ''],
    ['WheelLength', 'Angstrom', 1.0, [0, inf], '', ''],
    ['WheelSld', '', 1.0, [-inf, inf], '', ''],
    ['SolventSLD', '', 0.0, [-inf, inf], '', ''],
    ]


class GenSurface(object):
    def __init__(self, radius, diff):
        self.R = radius
        self.diff = diff

    umin = 0
    umax = np.pi
    vmin = 0
    vmax = 2 * np.pi

    def surface(self, u, v):
        """The default model is a sphere."""
        x = self.R * np.sin(u) * np.cos(v)
        y = self.R * np.sin(u) * np.sin(v)
        z = self.R * np.cos(u)
        return (x, y, z)

    def norm(self, u, v):
        """The default model is a sphere."""
        x = np.sin(u) * np.cos(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(u)
        return (x, y, z)

    def integrator(self, u, v):
        """The default model is a sphere."""
        return self.R**2 * np.sin(v)

    def calculate(self, v, u, k):
        # x, y, z = self.surface(u, v)
        # R = np.sqrt(x**2 + y**2 + z**2)
        # nx, ny, nz = self.norm(u, v)
        # dot = (x * nx + y * ny + z * nz)/R
        R = self.R
        dot = 1
        return self.diff * spherical_jn(1, k * R)/(k*R) * dot * self.integrator(u, v)

    def scatter(self, k):
        return dblquad(self.calculate, self.umin, self.umax,
                       lambda _: self.vmin, lambda _: self.vmax,
                       args=(k,))[0]


def Iq(x, AxleRadius, AxleLength, AxleSld, WheelRadius, WheelLength, WheelSld,
       SolventSLD):
    """Absolute scattering"""
    drho = (AxleSld - SolventSLD)

    surf = GenSurface(AxleRadius, drho)

    return surf.scatter(x)**2

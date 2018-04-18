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
from scipy.special import jn
from scipy.integrate import quad

name = "WheelAndAxle"
title = "User model for WheelAndAxle"
description = """"""

parameters = [
#   ["name", "units", default, [lower, upper], "type", "description"],
    ['AxleRadius', '', 1.0, [0, inf], '', ''],
    ['AxleLength', '', 1.0, [0, inf], '', ''],
    ['AxleSld', '', 1.0, [-inf, inf], '', ''],
    ['WheelRadius', '', 1.0, [0, inf], '', ''],
    ['WheelLength', '', 1.0, [0, inf], '', ''],
    ['WheelSld', '', 1.0, [-inf, inf], '', ''],
    ['SolventSLD', '', 1.0, [-inf, inf], '', ''],
    ]


def Iq(x, AxleRadius, AxleLength, AxleSld, WheelRadius, WheelLength, WheelSld,
       SolventSLD):
    """Absolute scattering"""
    drho = (AxleSld - SolventSLD)
    # result = 2 * np.pi * drho/x * quad(
    #     lambda z: jv(1, x * np.sqrt(np.sqrt(AxleRadius**2+z**2)*x)) *
    #     AxleRadius / np.sqrt(AxleRadius**2+z**2),
    #     0, AxleLength)

    wall = 2 * np.pi * drho / x * \
        quad(
            lambda z: jn(1, x * np.sqrt(AxleRadius**2+z**2)) *
            AxleRadius**2 / np.sqrt(AxleRadius**2+z**2),
            -AxleLength/2, AxleLength/2)[0]

    cap1 = 2 * np.pi * drho / x * \
        quad(
            lambda z: jn(1, x * np.sqrt((AxleLength/2)**2+z**2)) *
            z * AxleLength / np.sqrt((AxleLength/2)**2+z**2),
            0, AxleRadius)[0]

    V = np.pi * AxleRadius**2 * AxleLength

    return (wall + 2 * cap1)**2/V

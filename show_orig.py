# DOI: 10.1021/acs.jctc.2c01151

from math import cos, exp, pi
import numpy as np
from matplotlib import pyplot as plt

# (nm) cutoff distance
Rc = 1.35


def vdw(r, eps, r0, a, b):
    exponent = exp(f(r, r0, a, b) * (1 - r / r0))
    return eps * (exponent - 2 * exponent**0.5)


def f(r, r0, a, b):
    if r <= 0.5 * r0:
        return a * (3 + 2 * cos((2 * pi * r) / r0))
    elif r <= r0:
        return a
    else:
        return b


def main():
    data = {
        "WO-WO": [5.209, 0.594, 7.494, 9.355],
        "PEOT-PEOT": [4.000, 0.512, 10.288, 17.138],
        "PEOM-PEOT": [3.327, 0.512, 8.047, 17.138],
        "PEOM-PEOM": [2.767, 0.512, 6.294, 17.138],
        "PEGT-PEGT": [4.488, 0.521, 7.793, 13.169],
        "PEGT-PEOT": [2.881, 0.444, 8.954, 15.023],
        "PEGT-PEOM": [2.088, 0.517, 7.004, 9.149],
        "PEOT-WO": [4.739, 0.536, 8.781, 12.662],
        "PEOM-WO": [4.196, 0.536, 6.868, 12.662],
        "PEGT-WO": [4.508, 0.562, 7.642, 11.099],
    }
    for name, record in data.items():
        [eps, r0, a, b] = record
        r_list = np.arange(0.35, Rc, 0.01)
        vdw_list = [vdw(r, eps, r0, a, b) for r in r_list]

        plt.plot(r_list, vdw_list, label=name)

    plt.ylim(-8, 8)
    plt.legend(loc="upper right")
    plt.show()


main()

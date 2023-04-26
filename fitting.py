# DOI: 10.1021/acs.jctc.2c01151

from math import cos, exp, pi
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

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


def lj(r, eps, sigma):
    return 4 * eps * ((sigma / r) ** 12 - (sigma / r) ** 6)


def backingham(r, a, b, c):
    return a * np.exp(-b * r) - c / r**6


def main():
    plt.figure(figsize=(16, 9), dpi=120)

    data = {
        "WO-WO": [5.209, 0.594, 7.494, 9.355, "blue"],
        "PEOM-PEOM": [2.767, 0.512, 6.294, 17.138, "red"],
        "PEGT-PEGT": [4.488, 0.521, 7.793, 13.169, "green"],
        "PEGT-PEOM": [2.088, 0.517, 7.004, 9.149, "cyan"],
        "PEOM-WO": [4.196, 0.536, 6.868, 12.662, "magenta"],
        "PEGT-WO": [4.508, 0.562, 7.642, 11.099, "orange"],
    }
    for name, record in data.items():
        [eps, r0, a, b, color] = record
        r_list_morse = np.arange(0.3, Rc, 0.01)
        r_list_morse = list(filter(lambda x: vdw(x, eps, r0, a, b) < 8, r_list_morse))
        vdw_list_morse = [vdw(r, eps, r0, a, b) for r in r_list_morse]
        plt.plot(
            r_list_morse,
            vdw_list_morse,
            label=name + " (Morse)",
            color=color,
            linestyle="solid",
        )

        [eps, sigma], _ = curve_fit(lj, r_list_morse, vdw_list_morse)
        r_list_lj = np.copy(r_list_morse)
        vdw_list_lj = [lj(r, eps, sigma) for r in r_list_lj]
        plt.plot(
            r_list_lj,
            vdw_list_lj,
            label=name + " (LJ)",
            color=color,
            linestyle="dotted",
        )

        [a, b, c], _ = curve_fit(backingham, r_list_morse, vdw_list_morse)
        r_list_backingham = np.copy(r_list_morse)
        vdw_list_backingham = backingham(r_list_backingham, a, b, c)
        plt.plot(
            r_list_backingham,
            vdw_list_backingham,
            label=name + " (Backingham)",
            color=color,
            linestyle="dashed",
        )

    plt.ylim(-8, 8)
    plt.xlabel("$r_{ij}$ (nm)")
    plt.ylabel("$E_{pot}$ (kJ/mol)")
    plt.legend(bbox_to_anchor=(1, 1))
    plt.subplots_adjust(left=0.1, right=0.6, bottom=0.1, top=0.9)
    plt.savefig(
        "fitting.png",
    )


main()

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

plt.rcParams["font.size"] = "10.5"
plt.style.use("science")

naca_12 = [
    (0.0000000, 0.0000000),
    (0.0005839, 0.0042603),
    (0.0023342, 0.0084289),
    (0.0052468, 0.0125011),
    (0.0093149, 0.0164706),
    (0.0145291, 0.0203300),
    (0.0208771, 0.0240706),
    (0.0283441, 0.0276827),
    (0.0369127, 0.0311559),
    (0.0465628, 0.0344792),
    (0.0572720, 0.0376414),
    (0.0690152, 0.0406310),
    (0.0817649, 0.0434371),
    (0.0954915, 0.0460489),
    (0.1101628, 0.0484567),
    (0.1257446, 0.0506513),
    (0.1422005, 0.0526251),
    (0.1594921, 0.0543715),
    (0.1775789, 0.0558856),
    (0.1964187, 0.0571640),
    (0.2159676, 0.0582048),
    (0.2361799, 0.0590081),
    (0.2570083, 0.0595755),
    (0.2784042, 0.0599102),
    (0.3003177, 0.0600172),
    (0.3226976, 0.0599028),
    (0.3454915, 0.0595747),
    (0.3686463, 0.0590419),
    (0.3921079, 0.0583145),
    (0.4158215, 0.0574033),
    (0.4397317, 0.0563200),
    (0.4637826, 0.0550769),
    (0.4879181, 0.0536866),
    (0.5120819, 0.0521620),
    (0.5362174, 0.0505161),
    (0.5602683, 0.0487619),
    (0.5841786, 0.0469124),
    (0.6078921, 0.0449802),
    (0.6313537, 0.0429778),
    (0.6545085, 0.0409174),
    (0.6773025, 0.0388109),
    (0.6996823, 0.0366700),
    (0.7215958, 0.0345058),
    (0.7429917, 0.0323294),
    (0.7638202, 0.0301515),
    (0.7840324, 0.0279828),
    (0.8035813, 0.0258337),
    (0.8224211, 0.0237142),
    (0.8405079, 0.0216347),
    (0.8577995, 0.0196051),
    (0.8742554, 0.0176353),
    (0.8898372, 0.0157351),
    (0.9045085, 0.0139143),
    (0.9182351, 0.0121823),
    (0.9309849, 0.0105485),
    (0.9427280, 0.0090217),
    (0.9534372, 0.0076108),
    (0.9630873, 0.0063238),
    (0.9716559, 0.0051685),
    (0.9791229, 0.0041519),
    (0.9854709, 0.0032804),
    (0.9906850, 0.0025595),
    (0.9947532, 0.0019938),
    (0.9976658, 0.0015870),
    (0.9994161, 0.0013419),
    (1.0000000, 0.0012600),
]


# def coeff(n=4):

#     coeff[0] = 0.
#     # sum(coeff) = 0.
#     return coeff

coeff = np.random.rand(4)


def func(x, a, b, c, d, e, f, g, h, i, j):
    return (
        a * x
        + b * x**2
        + c * x**3
        + d * x**4
        + e * x**5
        + f * x**6
        + g * x**7
        + h * x**8
        + i * x**9
        + j * x**10
    )


def dfunc(x, a, b, c, d, e, f, g, h, i, j):
    return (
        a
        + 2 * b * x
        + 3 * c * x**2
        + 4 * d * x**3
        + 5 * e * x**4
        + 6 * f * x**5
        + 7 * g * x**6
        + 8 * h * x**7
        + 9 * i * x**8
        + 10 * j * x**9
    )


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def adjust(xdata):
    sc = 0.32
    adjust_factor = -sigmoid(3 + xdata * 15) * sc + sc
    return adjust_factor

def plot_data():
    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    ax.set_xlabel(r"x")
    ax.set_ylabel(r"y")

    ax.set_xlim(0, 1)
    # ax.set_ylim(0,0.1)

    xdata, ydata = np.array([p[0] for p in naca_12]), np.array([p[1] for p in naca_12])

    ax.plot(xdata, ydata, marker="+", label="NACA0012")

    # --- Adjust nose --- #
    ax.plot(xdata, adjust(xdata), label=r"$\sigma$ adjust factor")
    ydata = ydata - adjust(xdata)

    ax.plot(xdata, ydata, marker="*", label="Sharper nose")

    popt, pcov = curve_fit(func, xdata, ydata)
    print([round(pop, 3) for pop in popt])
    ax.plot(xdata, func(xdata, *popt), label=f"Polynomial $O({10})$")

    ax.legend()
    plt.savefig("foil_shape.pdf", dpi=200)


def plot_gradient():
    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    ax.set_xlabel(r"x")
    ax.set_ylabel(r"dy/dx")

    xdata, ydata = np.array([p[0] for p in naca_12]), np.array([p[1] for p in naca_12])

    ax.plot(xdata, np.gradient(ydata, xdata), label="NACA0012")
    ydata = ydata - adjust(xdata)

    popt, pcov = curve_fit(func, xdata, ydata)
    ax.plot(xdata, dfunc(xdata, *popt), label=f"Polynomial $O({10})$")

    ax.legend()
    plt.savefig("foil_grad.pdf", dpi=200)


if __name__ == "__main__":
    plot_data()
    plot_gradient()

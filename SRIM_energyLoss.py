import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from SRIM_Parser import SRIM_Table


# Constants
_c = 3.0e8              # m / s
_U = 1.66053906660E-24  # g
_avo = 6.02214076E23    # particles / mol

rho_24Mg = 1.74         # g / cm^2
mass_24Mg = 23.985041   # g / mol


def beamEnergyLoss(table,E):
    """ Returns the spline of the dE/dX for the beam species of the table """

    e_Knots = np.linspace(0,E,1000)
    dXdE_E = table.dXdE(e_Knots)
    dEdX_E = table.dEdX(e_Knots)

    # print(table.MaxE0)
    # print(table.MaxR)

    X = [0.]*len(e_Knots-1)
    step = len(e_Knots)-1-1

    # perform trapezoidal integration (b-a)*(f(b)+f(a))/2.
    while step >=0:
        trapezoid = X[step+1]+0.5*(e_Knots[step+1]-e_Knots[step]) * \
                    (dXdE_E[step+1]+dXdE_E[step])
        X[step] = trapezoid
        step -=1

    # Must reverse lists since CubicSpline requires ascending 'x' values
    X = np.asarray(list(reversed(X)))
    dEdX_E = np.asarray(list(reversed(dEdX_E)))
    spline = CubicSpline(X,dEdX_E)

    # convert units to MeV/um
    dEdX_E = dEdX_E/1000.

    # calculate polynomial
    z = np.polyfit(X,dEdX_E, 9)    # performs least squares polynomial fit of degree set
    f = np.poly1d(z)                # calculates the new points

    # print(z) # print the fit coefficients (order is highest order -> lowest order)

    # calculate new x's and y's
    x_new = np.linspace(X[0], X[-1], 1000)
    y_new = f(x_new)

    plt.scatter(x_new,y_new,s=1,alpha=.5,c='r')
    plt.plot(X,dEdX_E)
    plt.title("Bragg Curve")
    plt.xlabel("Range (um)")
    plt.ylabel("Stopping Power (MeV/um)")
    plt.savefig("Plots/He_in_Mg_braggCurve.png",dpi=900)
    # plt.show()
    plt.clf()


    # Convert to MeV
    e_Knots = e_Knots/1000

    # calculate polynomial
    zz = np.polyfit(e_Knots,dEdX_E, 9)  # performs least squares polynomial fit of degree set
    ff = np.poly1d(zz)                  # calculates the new points

    # print(z) # print the fit coefficients (order is highest order -> lowest order)

    # calculate new x's and y's
    xx_new = np.linspace(e_Knots[0], e_Knots[-1], 1000)
    yy_new = ff(xx_new)

    plt.scatter(xx_new,yy_new,s=1,alpha=.5,c='r')
    plt.plot(e_Knots,dEdX_E)
    plt.title("Stopping Crosssection")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Stopping Power (MeV/um)")
    plt.savefig("Plots/He_in_Mg_stoppingCross.png",dpi=900)

    return spline



"""
MAIN
"""
# Stopping Units = keV / micron
table_4He = SRIM_Table("SRIM_Data/Helium_in_Magnesium.txt")

ebeamMax = 5400
spline4He = beamEnergyLoss(table_4He,ebeamMax)  # up to 5.4 MeV beam energy converted to keV
eloss4He = spline4He.integrate(0,34.1)      # Caclulate energy loss due to length traveled by beam in um

# print(eloss4He)
# print(ebeamMax-eloss4He)

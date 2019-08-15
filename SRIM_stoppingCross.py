import numpy as np
import matplotlib.pyplot as plt
from SRIM_Parser import SRIM_Table


# Constants
_c = 3.0e8              # m / s
_U = 1.66053906660E-24  # g
_avo = 6.02214076E23    # particles / mol

# rho_24Mg = 1.74         # g / cm^2
# mass_24Mg = 23.985041   # g / mol


def stoppingCross(table,E,degree):
    """ Produces Stopping Crosssection: Stopping Power vs Energy """

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

    X = np.asarray(X)
    dEdX_E = np.asarray(dEdX_E)

    # convert units to MeV/um
    dEdX_E = dEdX_E/1000.

    # print(z) # print the fit coefficients (order is highest order -> lowest order)

    # Convert to MeV
    e_Knots = e_Knots/1000

    # calculate polynomial
    z = np.polyfit(e_Knots,dEdX_E, degree)  # performs least squares polynomial fit of degree set
    f = np.poly1d(z)                  # calculates the new points

    # print(z) # print the fit coefficients (order is highest order -> lowest order)

    # calculate new x's and y's
    x_new = np.linspace(e_Knots[0], e_Knots[-1], 1000)
    y_new = f(x_new)


    plt.plot(e_Knots,dEdX_E,c='k',label='Curve')
    plt.plot(x_new,y_new,alpha=1,c='r',label='Poly fit')
    plt.title("Stopping Crosssection")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Stopping Power (MeV/um)")
    plt.savefig("Plots/He_in_Mg_stoppingCross.png",dpi=900)


"""
MAIN
"""
# Stopping Units = keV / micron
table_4He = SRIM_Table("SRIM_Data/Helium_in_Magnesium.txt")

ebeamMax = 5400
stoppingCross(table_4He,ebeamMax,9)  # up to 5.4 MeV beam energy converted to keV

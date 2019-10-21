import numpy as np
import matplotlib.pyplot as plt
from SRIM_Parser import SRIM_Table


# Constants
_c = 3.0e8              # m / s
_U = 1.66053906660E-24  # g
_avo = 6.02214076E23    # particles / mol


def braggCurve(table,E,degree,color='k',label=''):
    """ Produces Bragg Curve: Stopping Power vs Range """

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
    # dEdX_E = dEdX_E/1000.
    dEdX_E = dEdX_E

    # calculate polynomial
    # z = np.polyfit(X,dEdX_E, degree)     # performs least squares polynomial fit of degree set
    # f = np.poly1d(z)                # calculates the new points

    # print(z) # print the fit coefficients (order is highest order -> lowest order)

    # calculate new x's and y's
    # x_new = np.linspace(X[0], X[-1], 1000)
    # y_new = f(x_new)

    # plt.plot(X,dEdX_E,c='k',label='Curve')
    plt.plot(X/1E4,dEdX_E*1E1,c=color,label=label)
    # plt.plot(x_new,y_new,alpha=1,c='r',label='Poly fit')
    plt.title("Bragg Curve")
    plt.xlabel("Range (cm)")
    # plt.ylabel("Stopping Power (MeV/um)")
    plt.ylabel("Stopping Power (MeV/cm)")
    # plt.legend()
    # plt.savefig("Plots/He_in_Mg_braggCurve.png",dpi=900)
    # plt.savefig("Plots/7Li_in_40Ar_braggCurve.png",dpi=900)
    # plt.savefig("Plots/8Li_in_40Ar_braggCurve.png",dpi=900)
    # plt.show()
    # plt.clf()

"""
MAIN
"""
# Stopping Units = keV / micron
# table_4He = SRIM_Table("SRIM_Data/Helium_in_Magnesium.txt")
table_40Ar_7 = SRIM_Table("SRIM_Data/7Lithium in Argon ~250Torr (gas).txt")
table_40Ar_8 = SRIM_Table("SRIM_Data/8Lithium in Argon ~250Torr (gas).txt")

# ebeamMax = 5400 # energy in keV
# braggCurve(table_4He,ebeamMax,9)  # up to 5.4 MeV beam energy converted to keV

ebeamMax = 21875
braggCurve(table_40Ar_7,ebeamMax,9,'r','$^7$Li')

ebeamMax = 25000
braggCurve(table_40Ar_8,ebeamMax,9,'b','$^8$Li')

plt.legend()
plt.savefig("Plots/Li_in_40Ar_braggCurves.png",dpi=900)

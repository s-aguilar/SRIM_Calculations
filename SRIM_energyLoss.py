import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from SRIM_Parser import SRIM_Table


# Constants
_c = 299792458          # m / s
_U = 1.66053906660E-24  # g
_avo = 6.02214076E23    # particles / mol

rho_24Mg = 1.74         # g / cm^3
mass_24Mg = 23.985041   # g / mol


def beamEnergyLoss(table,E,SU_conv):
    """ Returns the spline of the dE/dX for the beam species of the table,
        Energy 'E' is in keV, and 'SU' are the desired stopping units """

    e_Knots = np.linspace(0,E,1000)
    dXdE_E = table.dXdE(e_Knots)
    dEdX_E = table.dEdX(e_Knots)

    # print(table.MaxE0)
    # print(table.MaxR)

    # X = [0.]*len(e_Knots-1)

    # The index is 'step', remember arrays starts from zero
    # Integration of dX/dE_E vs E to obtain the X (distance of every energy step)
    # perform trapezoidal integration (b-a)*(f(b)+f(a))/2.
    X = np.zeros(len(e_Knots))
    step = len(e_Knots)-1-1
    while step >=0:
        trapezoid = X[step+1]+0.5*(e_Knots[step+1]-e_Knots[step]) * \
                    (dXdE_E[step+1]+dXdE_E[step])
        X[step] = trapezoid
        step -=1

    # Must reverse lists since CubicSpline requires ascending 'x' values
    X = np.asarray(list(reversed(X)))
    dEdX_E = np.asarray(list(reversed(dEdX_E)))
    spline = CubicSpline(X,dEdX_E)

    # convert units to desired stopping units
    dEdX_E = dEdX_E * SU_conv

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
    plt.xlabel("Range (%s)" % stopUnits.split()[2])
    plt.ylabel("Stopping Power (%s)" % stopUnits)
    # plt.tight_layout()
    plt.savefig("Plots/He_in_Mg_braggCurve.png",dpi=900)
    # plt.show()
    plt.clf()


    # Convert to Energy MeV
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
    plt.ylabel("Stopping Power (%s)" % stopUnits)
    # plt.tight_layout()
    plt.savefig("Plots/He_in_Mg_stoppingCross.png",dpi=900)

    return spline



"""
MAIN
"""

fName = 'SRIM_Data/Helium in Magnesium THICKNESS.txt'
table_4He = SRIM_Table(fName)
print("Stopping Units of SRIM file are:\n\t",table_4He.units)


# Choices for stopping Units:
#       'eV / A'
#       'keV / micron'
#       'MeV / mm'
#       'keV / (ug/cm2)'
#       'MeV / (mg/cm2)'
#       'keV / (mg/cm2)'
#       'eV / (1E15 atoms/cm2)'

stopUnits = 'keV / micron'
conv = table_4He.dictConversion[stopUnits]
print('\nDesired Stopping Units are:\n\t%s' % stopUnits)
print('\nA conversion factor of %f will be multiplied to the SRIM file\'s Stopping Units' % conv)
print('\n==================================================================\n')
ebeamMax = 4330

# Spline is dE/dX as function of X
spline4He = beamEnergyLoss(table_4He,ebeamMax,conv) # up to ebeamMax energy in keV, using conversion factor for desired units
eloss4He = spline4He.integrate(0,.14*conv)        # Caclulate energy loss due to length traveled by beam in um


print('For a %d keV beam, the energy loss is ' % ebeamMax, eloss4He,'keV')
# print(ebeamMax-eloss4He)

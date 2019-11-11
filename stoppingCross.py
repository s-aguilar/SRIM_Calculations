import numpy as np
import matplotlib.pyplot as plt
from SRIM_Parser import SRIM_Table

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


# Constants
_c = 299792458          # m / s
_U = 1.66053906660E-24  # g
_avo = 6.02214076E23    # particles / mol


def stoppingCross(table,E,low,high,SU_conv,stopUnitsDesired=False,eUnitsDesired=False,fit=False,degree=9,savePath='garbage.png'):
    """ Produces Stopping Crosssection: Stopping Power vs Beam Energy for the beam
    species of the table, Energy 'E' is in eV, 'low'and 'high' are define the beam
    energy range to fit over,'SU_conv' is a factor to multiply the default SRIM
    files stopping units into whatever 'stopUnitsDesired'. 'eUnitsDesired' is variable
    to select what energy units on the 'x' axis you want to use in the fit and the plot.
    'fit' is a boolean telling whether you would like to fit it with a 'degree'
    polynomial. 'savePath' is the location and name of which you would like to place
    the plot of the stopping cross section with or without the polynomial fit. """


    # Interpolate the stopping powers at these energies
    e_Knots = np.linspace(low,high,1000)

    # Stopping cross section curve points
    e_Knots_plot = np.linspace(0,ebeamMax,1000)

    # Convert stopping units if needed
    if SU_conv:
        dXdE_E = table.dXdE(e_Knots) * 1/SU_conv
        dEdX_E = table.dEdX(e_Knots) * SU_conv
        dXdE_E_plot = table.dXdE(e_Knots_plot) * 1/SU_conv
        dEdX_E_plot = table.dEdX(e_Knots_plot) * SU_conv
    else:
        dXdE_E = table.dXdE(e_Knots)
        dEdX_E = table.dEdX(e_Knots)
        dXdE_E_plot = table.dXdE(e_Knots_plot)
        dEdX_E_plot = table.dEdX(e_Knots_plot)

    if fit is True:

        # Calculate polynomial fit
        if eUnitsDesired:
            e_Knots *= eUnitsDesired[1]
            e_Knots_plot *= eUnitsDesired[1]
        z = np.polyfit(e_Knots,dEdX_E, degree)  # performs least squares polynomial fit of degree set
        f = np.poly1d(z)                        # calculates the new points

        print('\n\nPoly fit coefficients:\n[a%d, ..., a0]'%degree)
        print(z) # print the fit coefficients (order is highest order -> lowest order)

        # calculate stopping powers 'y'
        # x_new = np.linspace(e_Knots[0], e_Knots[-1], 1000)
        y_new = f(e_Knots)
        plt.plot(e_Knots,y_new,alpha=1,c='r',label='Poly fit order %d'%degree)
        plt.yscale('log')


    if eUnitsDesired:
        plt.xlabel('Energy ( %s )' % eUnitsDesired[0])
        plt.xlim(0,ebeamMax*eUnitsDesired[1])
    else:
        plt.xlabel('Energy ( %s )' % 'eV')
        plt.xlim(0,ebeamMax)

    plt.plot(e_Knots_plot,dEdX_E_plot,c='k',label='SRIM Stopping Curve',alpha=0.5)
    plt.title('Stopping Crosssection')
    plt.ylabel('Stopping Power ( %s )' % stopUnitsDesired)
    plt.legend()
    plt.savefig(savePath,dpi=900)
    plt.show()



"""
MAIN
"""

print('\n==================================================================')
print('==================================================================\n')

fName = 'SRIM_Data/Helium in Magnesium.txt'
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
#       'L.S.S. reduced units' NOT IMPLEMENTED

stopUnits = 'eV / (1E15 atoms/cm2)'
stopUnitsDesired = stopUnits


# Conversion factor multiplied to the default SRIM file's stopping units
# to attain 'stopUnits'
SU_conv = table_4He.dictConversion[stopUnits]


# If you want some other 'unorthodox' stopping units, manually adjust the
# conversion factor 'conv' and set unorthodox variable to 'True'
unorthodox = True


print('\nDesired SRIM Stopping Units are:\n\t%s' % stopUnits)

if not unorthodox:
    print('\nA conversion factor of %s will be multiplied to the SRIM \nfile\'s Stopping Units' % '{:.2e}'.format(SU_conv))
else:
    stopUnitsDesired = 'MeV / (atoms/cm2)'
    print('\nDesired Final Stopping Units are:\n\t%s' % stopUnitsDesired)

    # Unit conversion to energy you want
    # (SRIM stopping is in eV -> convert to MeV)
    E_conv = table_4He.dictOfUnits['MeV']

    # Unit conversion to Stopping you want
    U_conv = 1e15
    SU_conv = SU_conv * E_conv * U_conv

    print('\nA conversion factor of %s will instead be multiplied to the SRIM \nfile\'s Stopping Units' % '{:.2e}'.format(SU_conv))




print('\n==================================================================')
print('==================================================================\n')






# Max energy of your beam in eV
ebeamMax = 5.5e6

# 'low' energy of beam and 'high' in eV
low = 4e6
high = 5.45e6


print('Unit conversion factors used to attain desired units:\n',table_4He.dictOfUnits)


# 'False' if no fitting, 'True' if fitting with polynomial of order 'degree'
fit = True
degree = 3

# Desired Energy units for plotting
eUnitsDesired = ['MeV',table_4He.dictOfUnits['MeV']]
# eUnitsDesired = ['eV',table_4He.dictOfUnits['eV']]

savePath = 'Plots/He_in_Mg_stoppingCross.png'




print('\nPlotting ...')
stoppingCross(table_4He,ebeamMax,low,high,SU_conv,stopUnitsDesired,eUnitsDesired,fit,degree)
print('\nDone!\n\n')

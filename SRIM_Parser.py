import re
from scipy.interpolate import interp1d
import numpy as np

class SRIM_Table(object):
    def __init__(self,fname):

        # TO DO: 1. Speed up this process?
        #        2. Handle edge cases

        sci = '(\d+\.\d+E[+-]*\d+)' # check sign of power
        flo = '(\d+\.\d+|\d+)'      # '|' is OR. Check whether it is a float or an int
        unit = '(\w+)'
        fmt = ('\s*'+flo+'\s'+unit+'\s+'+sci+'\s+'+sci
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'.*')

        data = []

        # eV / A
        units_eVa = {'eV':1, 'keV':1e3, 'MeV':1e6, \
                        'm':1e10, 'mm':1e7, 'um':1e4, 'A':1}

        # keV / um
        units_keVum = {'eV':1e-3, 'keV':1, 'MeV':1e3, \
                        'm':1e-6, 'mm':1e-3, 'um':1, 'A':1e-4}

        # MeV / mm
        units_MeVmm = {'eV':1e-6, 'keV':1e-3, 'MeV':1, \
                        'm':1e3, 'mm':1, 'um':1e-3, 'A':1e-7}

        # keV / (ug/cm2)
        units_keVugcm2 = {'eV':1e-3, 'keV':1, 'MeV':1e3, \
                        'm':1e-2, 'mm':1e1, 'um':1e4, 'A':1e8}

        # MeV / (mg/cm2)
        units_MeVmgcm2 = {'eV':1e-6, 'keV':1e-3, 'MeV':1, \
                        'm':1e-2, 'mm':1e1, 'um':1e4, 'A':1e8}

        # keV / (mg/cm2)
        units_keVmgcm2 = {'eV':1e-3, 'keV':1, 'MeV':1e3, \
                        'm':1e-2, 'mm':1e1, 'um':1e4, 'A':1e8}

        # eV / (1E15 atoms/cm^2)
        units_eVatomscm2 = {'eV':1, 'keV':1e3, 'MeV':1e6, \
                        'm':1e-2, 'mm':1e1, 'um':1e4, 'A':1e8}


        # dictionary of units
        dictUnits = {'eV / A':units_eVa,'keV / micron':units_keVum,
            'MeV / mm':units_MeVmm,'keV / (ug/cm2)':units_keVugcm2,
            'MeV / (mg/cm2)':units_MeVmgcm2,'keV / (mg/cm2)':units_keVmgcm2,
            'eV / (1E15 atoms/cm2)':units_eVatomscm2}


        with open(fname) as f:

            def convert(line):
                """
                Clean up the line and extract conversion factor
                """

                lineSplit = line.split()
                return np.float64(lineSplit[0])


            for line in f:
                res = re.match(fmt,line)

                # Acquire the stopping units
                stopUnit = re.search(' Stopping Units = ',line)
                if stopUnit:
                    units = str(line[19:-1]).strip()    # strip white space
                    dictOfUnits = dictUnits[units]
                    continue

                # Acquire the conversion table between units
                conv = re.search(' Multiply Stopping by        for Stopping Units',line)
                if conv:
                    f.readline()    # Read empty line -------------------        ------------------
                    dictConversion = {'eV / A':convert(f.readline()),
                        'keV / micron':convert(f.readline()),
                        'MeV / mm':convert(f.readline()),
                        'keV / (ug/cm2)':convert(f.readline()),
                        'MeV / (mg/cm2)':convert(f.readline()),
                        'keV / (mg/cm2)':convert(f.readline()),
                        'eV / (1E15 atoms/cm2)':convert(f.readline()) }
                    # print(f.readline())   # L.S.S reduced units not interested in these
                    break

                elif res is None: continue

                # The '-u' suffix denotes unit and L is length/range
                E,Eu, dE_dX1, dE_dX2, L, Lu = res.groups()[:6]

                # print(E,Eu,dE_dX1, dE_dX2, L, Lu)
                Total_dEdX = np.float64(dE_dX1) + np.float64(dE_dX2)
                data.append({'E':np.float64(E)*dictOfUnits[Eu],
                'Range':np.float64(L)*dictOfUnits[Lu],
                'dXdE':1./Total_dEdX,'dEdX':Total_dEdX})

        if not data:
            raise Exception("Input a valid table")

        self.MaxE0 = max(_['E'] for _ in data)
        self.MaxR = max(_['Range'] for _ in data)
        self.dXdE = interp1d(*zip(*[(_['E'],_['dXdE'],) for _ in data]),
                fill_value='extrapolate' )
        self.dEdX = interp1d(*zip(*[(_['E'],_['dEdX'],) for _ in data]),
                fill_value='extrapolate' )
        self.R2E = interp1d(*zip(*[(_['Range'],_['E'],) for _ in data]),
                fill_value='extrapolate' )
        self.E2R = interp1d(*zip(*[(_['E'],_['Range'],) for _ in data]),
                fill_value='extrapolate' )
        self.units = units
        self.dictOfUnits = dictOfUnits
        self.dictConversion = dictConversion


    def GetE(self,E0,X):
        if np.logical_or.reduce(np.array(E0)>self.MaxE0):
            raise Exception("exceed the maximum energy range of the table")

        R = self.E2R(E0)
        R_ = R-X
        E = self.R2E(R_)
        return np.clip(E,0,float('inf'))

    def GetE0(self,x):
        if np.logical_or.reduce(np.array(x)>self.MaxR):
            raise Exception("exceed the maximum distance range of the table")

        return self.R2E(x)

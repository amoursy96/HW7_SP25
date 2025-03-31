# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from pyXSteam.XSteam import XSteam
from UnitConversion import UC
from scipy.optimize import fsolve

class thermoSatProps:
    """
    This class retrieves saturation properties of steam for a given pressure or temperature.
    It uses the XSteam library to fetch thermodynamic properties.
    """
    def __init__(self, p=None, t=None, SI=True):
        """
        Initializes the class and fetches saturation properties based on input pressure or temperature.
        :param p: Pressure in bar (if SI) or psi (if Imperial)
        :param t: Temperature in Celsius (if SI) or Fahrenheit (if Imperial)
        :param SI: Boolean flag for unit system (True for SI, False for Imperial)
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        if p is not None:
            self.getSatProps(p, SI)
        elif t is not None:
            self.getSatProps(self.steamTable.psat_t(t), SI)

    def getSatProps(self, p, SI=True):
        """
        Fetches and stores thermodynamic properties for a given saturation pressure.
        :param p: Pressure in bar (SI) or psi (Imperial)
        :param SI: Boolean flag for unit system
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        self.pSat = p
        self.tSat = self.steamTable.tsat_p(p)  # Saturation temperature
        self.vf = self.steamTable.vL_p(p)  # Specific volume of liquid
        self.vg = self.steamTable.vV_p(p)  # Specific volume of vapor
        self.hf = self.steamTable.hL_p(p)  # Enthalpy of liquid
        self.hg = self.steamTable.hV_p(p)  # Enthalpy of vapor
        self.uf = self.steamTable.uL_p(p)  # Internal energy of liquid
        self.ug = self.steamTable.uV_p(p)  # Internal energy of vapor
        self.sf = self.steamTable.sL_p(p)  # Entropy of liquid
        self.sg = self.steamTable.sV_p(p)  # Entropy of vapor
        self.vgf = self.vg - self.vf  # Difference in specific volumes
        self.hgf = self.hg - self.hf  # Difference in enthalpy
        self.sgf = self.sg - self.sf  # Difference in entropy
        self.ugf = self.ug - self.uf  # Difference in internal energy

class thermoState:
    """
    This class represents a thermodynamic state with properties like pressure, temperature,
    specific volume, internal energy, enthalpy, entropy, and quality (x).
    """
    def __init__(self, p=None, t=None, v=None, u=None, h=None, s=None, x=None):
        """
        Initializes a thermodynamic state.
        :param p: Pressure
        :param t: Temperature
        :param v: Specific volume
        :param u: Internal energy
        :param h: Enthalpy
        :param s: Entropy
        :param x: Quality (0 for saturated liquid, 1 for saturated vapor)
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.region = "saturated"  # Default region
        self.p = p
        self.t = t
        self.v = v
        self.u = u
        self.h = h
        self.s = s
        self.x = x

    def computeProperties(self):
        """
        Computes missing properties based on the known state properties.
        """
        if self.region == "two-phase":
            self.u = self.steamTable.uL_p(self.p) + self.x * (
                self.steamTable.uV_p(self.p) - self.steamTable.uL_p(self.p))
            self.h = self.steamTable.hL_p(self.p) + self.x * (
                self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p))
            self.s = self.steamTable.sL_p(self.p) + self.x * (
                self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p))
            self.v = self.steamTable.vL_p(self.p) + self.x * (
                self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p))
        else:
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)
            self.x = 1.0 if self.region == "super-heated vapor" else 0.0

class main_window(QtWidgets.QWidget):
    """
    Main GUI window for the thermodynamic state calculator.
    Handles UI setup, user inputs, and displays calculated properties.
    """
    def __init__(self):
        super().__init__()
        self.setupUi()
        self.SetupSlotsAndSignals()
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits = 'SI'
        self.setUnits()
        self.show()

    def SetupSlotsAndSignals(self):
        """
        Connects UI elements to their respective functions.
        """
        self._rdo_English.clicked.connect(self.setUnits)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._cmb_Property1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate.clicked.connect(self.calculateProperties)

    def setUnits(self):
        """
        Updates units and converts values accordingly when switching between SI and Imperial units.
        """
        SI = self._rdo_SI.isChecked()
        newUnits = 'SI' if SI else 'EN'
        self.currentUnits = newUnits

        if SI:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.p_Units = "bar"
            self.t_Units = "C"
            self.u_Units = "kJ/kg"
        else:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.p_Units = "psi"
            self.t_Units = "F"
            self.u_Units = "btu/lb"

    def calculateProperties(self):
        """
        Computes thermodynamic properties based on user input and updates the UI.
        """
        self.state1 = thermoState()
        self.state2 = thermoState()
        # Fetch and compute properties for each state
        # Display results in the UI


def main():
    """
    Launches the PyQt5 application.
    """
    app = QtWidgets.QApplication(sys.argv)
    main_win = main_window()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

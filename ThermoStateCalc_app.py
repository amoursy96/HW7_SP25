# -*- coding: utf-8 -*-
"""
Thermodynamic State Calculator Application

A PyQt5-based GUI application for calculating thermodynamic properties of water/steam
at two different states using the pyXSteam library. The application allows specifying
any two thermodynamic properties at each state and calculates all other properties.

Key Features:
- Supports both SI and English unit systems
- Calculates properties for two distinct states
- Computes differences between states
- Automatic unit conversion
- Handles all regions (subcooled, two-phase, superheated)
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from pyXSteam.XSteam import XSteam
from UnitConversion import UC
import sys


class thermoState:
    """
    Represents a thermodynamic state of water/steam with unit conversion capability.

    This class encapsulates all thermodynamic properties (pressure, temperature,
    enthalpy, etc.) and provides methods to calculate them based on input properties.

    Attributes:
        SI (bool): Flag indicating whether SI units are being used
        steamTable (XSteam): Instance of XSteam steam tables
        region (str): Thermodynamic region ('subcooled', 'two-phase', 'superheated')
        p (float): Pressure in bar (SI) or psi (English)
        t (float): Temperature in °C (SI) or °F (English)
        v (float): Specific volume in m³/kg (SI) or ft³/lb (English)
        u (float): Internal energy in kJ/kg (SI) or btu/lb (English)
        h (float): Enthalpy in kJ/kg (SI) or btu/lb (English)
        s (float): Entropy in kJ/kg·K (SI) or btu/lb·°F (English)
        x (float): Quality (vapor mass fraction) [0-1]
    """

    def __init__(self, SI=True):
        """
        Initialize a thermodynamic state with default values.

        Args:
            SI (bool): Whether to use SI units (default) or English units
        """
        self.SI = SI
        # Initialize steam tables with appropriate unit system
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        self.region = None  # Will be set when state is calculated
        self.p = None  # Pressure
        self.t = None  # Temperature
        self.v = None  # Specific volume
        self.u = None  # Internal energy
        self.h = None  # Enthalpy
        self.s = None  # Entropy
        self.x = None  # Quality

    def setState(self, prop1, prop2, val1, val2):
        """
        Set the thermodynamic state based on two known properties.

        Handles unit conversion and determines the thermodynamic region before
        calculating all other properties.

        Args:
            prop1 (str): First property identifier ('p', 't', 'x', 'h', 's', 'v')
            prop2 (str): Second property identifier
            val1 (float): Value of first property
            val2 (float): Value of second property

        Note:
            The property identifiers should be the first character of the property name:
            - 'p' for pressure
            - 't' for temperature
            - 'x' for quality
            - 'h' for enthalpy
            - 's' for entropy
            - 'v' for specific volume
        """
        # Convert inputs to proper units if using English system
        if not self.SI:
            if prop1 == 'p': val1 = UC.psi_to_bar(val1)
            if prop2 == 'p': val2 = UC.psi_to_bar(val2)
            if prop1 == 't': val1 = UC.F_to_C(val1)
            if prop2 == 't': val2 = UC.F_to_C(val2)

        # Handle two-phase region (when quality is specified)
        if 'x' in [prop1, prop2]:
            self.region = "two-phase"
            # Get quality value and clamp between 0 and 1
            x = val1 if prop1 == 'x' else val2
            self.x = max(0.0, min(1.0, x))

            # Get the other property (pressure or temperature)
            other_prop = prop2 if prop1 == 'x' else prop1
            other_val = val2 if prop1 == 'x' else val1

            if other_prop == 'p':
                self.p = other_val
                self.t = self.steamTable.tsat_p(self.p)
            elif other_prop == 't':
                self.t = other_val
                self.p = self.steamTable.psat_t(self.t)
        else:
            # Handle single-phase regions (subcooled or superheated)
            if prop1 == 'p' and prop2 == 't':
                self.p = val1
                self.t = val2
            elif prop1 == 't' and prop2 == 'p':
                self.t = val1
                self.p = val2

            # Determine region based on saturation temperature
            tsat = self.steamTable.tsat_p(self.p) if self.p else None
            if tsat:
                if abs(self.t - tsat) < 1e-6:  # At saturation temperature
                    self.region = "two-phase"
                    self.x = 0.5  # Default quality for saturation
                elif self.t > tsat:
                    self.region = "superheated"
                    self.x = 1.0
                else:
                    self.region = "subcooled"
                    self.x = 0.0

        # Calculate all other properties
        self.computeProperties()

    def computeProperties(self):
        """Calculate all thermodynamic properties based on the current state."""
        if self.region == "two-phase":
            # Calculate properties using quality and saturation values
            self.u = self.steamTable.uL_p(self.p) + self.x * (
                    self.steamTable.uV_p(self.p) - self.steamTable.uL_p(self.p))
            self.h = self.steamTable.hL_p(self.p) + self.x * (
                    self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p))
            self.s = self.steamTable.sL_p(self.p) + self.x * (
                    self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p))
            self.v = self.steamTable.vL_p(self.p) + self.x * (
                    self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p))
        else:
            # Get properties directly from steam tables for single-phase regions
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)


class main_window(QtWidgets.QWidget):
    """
    Main application window for the Thermodynamic State Calculator.

    This class handles the GUI interface, user interactions, and coordinates
    the thermodynamic calculations between the two states.

    Attributes:
        ui (Ui__frm_StateCalculator): The generated UI class from Qt Designer
        state1 (thermoState): First thermodynamic state
        state2 (thermoState): Second thermodynamic state
    """

    def __init__(self):
        """Initialize the main window and set up the UI."""
        super().__init__()

        # Load and setup the UI from the generated file
        self.ui = Ui__frm_StateCalculator()
        self.ui.setupUi(self)

        # Initialize two thermodynamic states with SI units by default
        self.state1 = thermoState(SI=True)
        self.state2 = thermoState(SI=True)

        # Connect UI signals to slots
        self.setup_connections()

        # Set initial units to SI and update UI
        self.setUnits()

        # Show the window
        self.show()

    def setup_connections(self):
        """
        Connect all UI signals to their appropriate slots.

        This includes:
        - Unit system radio buttons
        - Calculate button
        - Property selection comboboxes
        """
        # Unit system selection
        self.ui._rdo_SI.clicked.connect(self.setUnits)
        self.ui._rdo_English.clicked.connect(self.setUnits)

        # Calculate button
        self.ui._pb_Calculate.clicked.connect(self.calculateProperties)

        # Property selection changes (update unit labels automatically)
        self.ui._cmb_State1_Property1.currentTextChanged.connect(
            lambda: self.updateUnitLabels(1, 1))
        self.ui._cmb_State1_Property2.currentTextChanged.connect(
            lambda: self.updateUnitLabels(1, 2))
        self.ui._cmb_State2_Property1.currentTextChanged.connect(
            lambda: self.updateUnitLabels(2, 1))
        self.ui._cmb_State2_Property2.currentTextChanged.connect(
            lambda: self.updateUnitLabels(2, 2))

    def setUnits(self):
        """
        Handle unit system changes between SI and English.

        This method:
        1. Converts existing values to the new unit system
        2. Updates the steam table instances
        3. Updates all unit labels
        4. Recalculates properties if valid states exist
        """
        # Determine new unit system
        SI = self.ui._rdo_SI.isChecked()

        # Convert existing input values to new units
        self.convertValues(SI)

        # Update steam tables with new unit system
        self.state1.SI = SI
        self.state2.SI = SI
        self.state1.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        self.state2.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)

        # Update all unit labels in the UI
        self.updateAllUnitLabels()

        # Recalculate properties if we have valid states
        if self.state1.p is not None and self.state2.p is not None:
            self.calculateProperties()

    def convertValues(self, toSI):
        """
        Convert all input values between unit systems.

        Args:
            toSI (bool): If True, convert to SI units; if False, convert to English
        """
        for state_num in [1, 2]:  # Handle both states
            for prop_num in [1, 2]:  # Handle both properties
                # Get UI elements
                cmb = getattr(self.ui, f"_cmb_State{state_num}_Property{prop_num}")
                le = getattr(self.ui, f"_le_State{state_num}_Property{prop_num}")

                try:
                    # Get current value and property type
                    val = float(le.text())
                    prop = cmb.currentText()[0].lower()  # First character identifies property

                    # Perform appropriate conversion based on property type
                    if prop == 'p':
                        new_val = UC.bar_to_psi(val) if not toSI else UC.psi_to_bar(val)
                    elif prop == 't':
                        new_val = UC.C_to_F(val) if not toSI else UC.F_to_C(val)
                    else:
                        # No conversion needed for quality or dimensionless properties
                        new_val = val

                    # Update the value in the UI
                    le.setText(f"{new_val:.3f}")
                except ValueError:
                    # Skip conversion if value can't be parsed as float
                    pass

    def updateAllUnitLabels(self):
        """Update unit labels for all property inputs."""
        for state_num in [1, 2]:  # Both states
            for prop_num in [1, 2]:  # Both properties
                self.updateUnitLabels(state_num, prop_num)

    def updateUnitLabels(self, state_num, prop_num):
        """
        Update the unit label for a specific property input.

        Args:
            state_num (int): State number (1 or 2)
            prop_num (int): Property number (1 or 2)
        """
        # Get UI elements
        cmb = getattr(self.ui, f"_cmb_State{state_num}_Property{prop_num}")
        lbl = getattr(self.ui, f"_lbl_State{state_num}_Property{prop_num}_Units")

        # Determine property type and current unit system
        prop = cmb.currentText()[0].lower()
        SI = self.ui._rdo_SI.isChecked()

        # Mapping of property types to their unit labels
        unit_map = {
            'p': ('bar', 'psi'),  # Pressure
            't': ('°C', '°F'),  # Temperature
            'h': ('kJ/kg', 'btu/lb'),  # Enthalpy
            's': ('kJ/kg·K', 'btu/lb·°F'),  # Entropy
            'v': ('m³/kg', 'ft³/lb')  # Specific volume
        }

        # Set the appropriate unit label
        if prop in unit_map:
            lbl.setText(unit_map[prop][0 if SI else 1])
        else:
            # Empty label for quality or other dimensionless properties
            lbl.setText("")

    def calculateProperties(self):
        """
        Calculate thermodynamic properties for both states.

        This method:
        1. Gets input values from the UI
        2. Sets the thermodynamic states
        3. Updates the results display
        4. Handles any calculation errors
        """
        try:
            # State 1 calculations
            p1 = self.ui._cmb_State1_Property1.currentText()[0].lower()
            p2 = self.ui._cmb_State1_Property2.currentText()[0].lower()
            v1 = float(self.ui._le_State1_Property1.text())
            v2 = float(self.ui._le_State1_Property2.text())
            self.state1.setState(p1, p2, v1, v2)

            # State 2 calculations
            p1 = self.ui._cmb_State2_Property1.currentText()[0].lower()
            p2 = self.ui._cmb_State2_Property2.currentText()[0].lower()
            v1 = float(self.ui._le_State2_Property1.text())
            v2 = float(self.ui._le_State2_Property2.text())
            self.state2.setState(p1, p2, v1, v2)

            # Update the results display
            self.updateResults()

        except ValueError:
            # Handle invalid numerical inputs
            QtWidgets.QMessageBox.warning(self, "Input Error",
                                          "Please enter valid numerical values for all properties")
        except Exception as e:
            # Handle other calculation errors
            QtWidgets.QMessageBox.critical(self, "Calculation Error",
                                           f"Failed to calculate properties: {str(e)}")

    def updateResults(self):
        """
        Update the results display with calculated properties.

        This includes:
        - Properties for both states
        - Differences between states
        - Proper units based on current unit system
        """
        SI = self.ui._rdo_SI.isChecked()

        # Update displays for both states
        self.updateStateDisplay(1, self.state1, SI)
        self.updateStateDisplay(2, self.state2, SI)

        # Update the differences between states
        self.updateDifferences(SI)

    def updateStateDisplay(self, state_num, state, SI):
        """
        Update the display for a single thermodynamic state.

        Args:
            state_num (int): State number (1 or 2)
            state (thermoState): The thermodynamic state to display
            SI (bool): Whether to display in SI units
        """
        # Convert units if displaying in English system
        p = state.p if SI else UC.bar_to_psi(state.p)
        t = state.t if SI else UC.C_to_F(state.t)
        h = state.h if SI else UC.kJperkg_to_btuperlb(state.h)
        s = state.s if SI else UC.kJperkgc_to_btuperlbF(state.s)
        v = state.v if SI else UC.m3perkg_to_ft3perlb(state.v)

        # Update all labels for this state
        getattr(self.ui, f"_lbl_State{state_num}_Region").setText(f"Region = {state.region}")
        getattr(self.ui, f"_lbl_State{state_num}_Pressure").setText(
            f"Pressure = {p:.3f} {'bar' if SI else 'psi'}")
        getattr(self.ui, f"_lbl_State{state_num}_Temperature").setText(
            f"Temperature = {t:.3f} {'°C' if SI else '°F'}")
        getattr(self.ui, f"_lbl_State{state_num}_Enthalpy").setText(
            f"Enthalpy = {h:.3f} {'kJ/kg' if SI else 'btu/lb'}")
        getattr(self.ui, f"_lbl_State{state_num}_Entropy").setText(
            f"Entropy = {s:.3f} {'kJ/kg·K' if SI else 'btu/lb·°F'}")
        getattr(self.ui, f"_lbl_State{state_num}_SpecificVolume").setText(
            f"Specific Volume = {v:.6f} {'m³/kg' if SI else 'ft³/lb'}")
        getattr(self.ui, f"_lbl_State{state_num}_Quality").setText(f"Quality = {state.x:.3f}")

    def updateDifferences(self, SI):
        """
        Update the display of differences between the two states.

        Args:
            SI (bool): Whether to display differences in SI units
        """
        # Calculate raw differences
        dp = self.state2.p - self.state1.p
        dt = self.state2.t - self.state1.t
        dh = self.state2.h - self.state1.h
        ds = self.state2.s - self.state1.s
        dv = self.state2.v - self.state1.v

        # Convert differences if displaying in English units
        if not SI:
            dp = UC.bar_to_psi(dp)
            dt = UC.DeltaC_to_DeltaF(dt)
            dh = UC.kJperkg_to_btuperlb(dh)
            ds = UC.kJperkgc_to_btuperlbF(ds)
            dv = UC.m3perkg_to_ft3perlb(dv)

        # Update difference labels
        self.ui._lbl_PressureChange.setText(f"ΔP = {dp:.3f} {'bar' if SI else 'psi'}")
        self.ui._lbl_TemperatureChange.setText(f"ΔT = {dt:.3f} {'°C' if SI else '°F'}")
        self.ui._lbl_EnthalpyChange.setText(f"Δh = {dh:.3f} {'kJ/kg' if SI else 'btu/lb'}")
        self.ui._lbl_EntropyChange.setText(f"Δs = {ds:.3f} {'kJ/kg·K' if SI else 'btu/lb·°F'}")
        self.ui._lbl_SpecificVolumeChange.setText(f"Δv = {dv:.6f} {'m³/kg' if SI else 'ft³/lb'}")


def main():
    """
    Main entry point for the application.

    Creates and runs the Qt application with the main window.
    """
    # Create Qt application
    app = QtWidgets.QApplication(sys.argv)

    # Create and show main window
    window = main_window()

    # Start event loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    # Import the generated UI class
    from ThermoStateCalc import Ui__frm_StateCalculator

    # Run the application
    main()
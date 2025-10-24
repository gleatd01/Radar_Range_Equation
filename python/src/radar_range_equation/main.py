import sympy
import sympy.physics.units as units
from sympy.physics.units import convert_to
from sympy.physics.units import speed_of_light

class Variables:
    """
    A class to hold variables that can be dynamically assigned and redefined.
    """
    def __init__(self, **overrides):
        # Default physical constants and other useful defaults
        self.c = convert_to(speed_of_light, [units.meter, units.seconds])          # speed of light (m/s)
        self.k = 1.38e-23           # Boltzmann constant (J/K) — example
        self.T_0 = 290               # reference temperature (K) — example
        # prefer non-keyword name for wavelength; keep for convenience too
        self.lambda_ = None         # preferred attribute name for wavelength
        # apply any overrides passed in
        for name, val in overrides.items():
            setattr(self, name, val)
    pass

vars = Variables()


def redefine_variable(var_name, new_value):
    """
    Redefines a global variable within the 'vars' namespace.
    Args:
        var_name (str): The name of the variable to redefine (e.g., "lambda").
        new_value: The new value to assign to the variable.
    """
    setattr(vars, var_name, new_value)

# Demonstration of usage from a separate script or module:
if __name__ == '__main__':  # Only runs when the script is executed directly
    import radar_range_equation as RRE
    # Use setattr to set 'lambda' since it's a reserved keyword
    setattr(RRE.vars, 'lambda', .03)
    
    print(f"Initial lambda (RRE.vars.lambda): {getattr(RRE.vars, 'lambda')}")

    RRE.redefine_variable("lambda", 0.07)  # Redefine it through the package

    print(f"Redefined lambda (RRE.vars.lambda): {getattr(RRE.vars, 'lambda')}")

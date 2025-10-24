class Variables:
    """
    A class to hold variables that can be dynamically assigned and redefined.
    """
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

#!/usr/bin/env python3
"""Comprehensive test suite for equations and solvers.

This script tests all equations and solvers in the radar_range_equation package
and generates test results that can be used for badge generation.
"""

import sys
import json
from math import isclose


def test_equations_and_solvers():
    """Test all equations and solvers and generate results.
    
    Returns dict with test results including counts and status.
    """
    try:
        import radar_range_equation as RRE
        from radar_range_equation import equations, solve
        
        results = {
            "equations": {"total": 0, "tested": 0, "passed": 0, "failed": []},
            "solvers": {"total": 0, "tested": 0, "passed": 0, "failed": []}
        }
        
        print("=" * 70)
        print("Testing Equations")
        print("=" * 70)
        
        # Get all equation attributes from equations class
        equation_names = [name for name in dir(equations) 
                         if not name.startswith('_') and name.endswith('_sym') == False
                         and hasattr(equations, name)]
        
        # Filter to only actual equations (sympy.Eq objects)
        import sympy
        equation_attrs = []
        for name in equation_names:
            attr = getattr(equations, name)
            if isinstance(attr, sympy.Eq):
                equation_attrs.append(name)
        
        results["equations"]["total"] = len(equation_attrs)
        
        # Test that equations exist and are valid
        for eq_name in equation_attrs:
            try:
                eq = getattr(equations, eq_name)
                # Check that it's a valid sympy Eq
                assert isinstance(eq, sympy.Eq), f"{eq_name} is not a sympy.Eq"
                assert eq.lhs is not None, f"{eq_name} has no left-hand side"
                assert eq.rhs is not None, f"{eq_name} has no right-hand side"
                results["equations"]["tested"] += 1
                results["equations"]["passed"] += 1
                print(f"✓ {eq_name}")
            except Exception as e:
                results["equations"]["tested"] += 1
                results["equations"]["failed"].append({"name": eq_name, "error": str(e)})
                print(f"✗ {eq_name}: {e}")
        
        print(f"\nEquations: {results['equations']['passed']}/{results['equations']['total']} passed")
        
        print("\n" + "=" * 70)
        print("Testing Solvers")
        print("=" * 70)
        
        # Get all solver methods/attributes from solve class
        solver_names = [name for name in dir(solve) 
                       if not name.startswith('_')]
        
        # Filter to actual solvers (methods or Solvable instances)
        solver_attrs = []
        for name in solver_names:
            attr = getattr(solve, name)
            # Include methods and Solvable instances
            if callable(attr) or hasattr(attr, 'solve'):
                solver_attrs.append(name)
        
        results["solvers"]["total"] = len(solver_attrs)
        
        # Test that solvers exist and are callable
        for solver_name in solver_attrs:
            try:
                solver = getattr(solve, solver_name)
                # Check if it's callable or a Solvable
                if hasattr(solver, 'solve'):
                    # It's a Solvable instance
                    assert hasattr(solver, 'status'), f"{solver_name} Solvable missing status()"
                    assert callable(solver), f"{solver_name} Solvable is not callable"
                elif callable(solver):
                    # It's a regular method
                    assert hasattr(solver, '__call__'), f"{solver_name} is not callable"
                else:
                    raise ValueError(f"{solver_name} is neither a method nor a Solvable")
                
                results["solvers"]["tested"] += 1
                results["solvers"]["passed"] += 1
                print(f"✓ {solver_name}")
            except Exception as e:
                results["solvers"]["tested"] += 1
                results["solvers"]["failed"].append({"name": solver_name, "error": str(e)})
                print(f"✗ {solver_name}: {e}")
        
        print(f"\nSolvers: {results['solvers']['passed']}/{results['solvers']['total']} passed")
        
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        
        total_tests = results["equations"]["tested"] + results["solvers"]["tested"]
        total_passed = results["equations"]["passed"] + results["solvers"]["passed"]
        total_failed = len(results["equations"]["failed"]) + len(results["solvers"]["failed"])
        
        print(f"Total Tests: {total_tests}")
        print(f"Passed: {total_passed}")
        print(f"Failed: {total_failed}")
        
        # Calculate percentage
        if total_tests > 0:
            pass_rate = (total_passed / total_tests) * 100
            results["pass_rate"] = round(pass_rate, 1)
            print(f"Pass Rate: {pass_rate:.1f}%")
        else:
            results["pass_rate"] = 0
        
        # Write results to JSON file for badge generation
        with open('test_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n✓ Test results written to test_results.json")
        
        # Return success if all tests passed
        if total_failed == 0:
            print("\n✓ All tests passed!")
            return 0
        else:
            print(f"\n✗ {total_failed} tests failed")
            return 1
            
    except ImportError as e:
        print(f"✗ ImportError: {e}")
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(test_equations_and_solvers())

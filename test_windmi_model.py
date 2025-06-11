"""
Unit tests for the WINDMI model implementation in `windmi_model.py`.

This test suite uses the `unittest` framework to verify the correctness
of the `solve_windmi_model` function, including its behavior with
example parameters, specific scenarios like zero input, and the
expected shape and type of its output.
"""
import unittest
import numpy as np
from windmi_model import solve_windmi_model, example_p_in, example_p_out as example_p_out_original

class TestWindmiModel(unittest.TestCase):
    """
    Test class for the WINDMI model solver.

    Contains various test methods to ensure the reliability of the
    `solve_windmi_model` function under different conditions.
    """

    def setUp(self):
        """Set up common parameters for tests."""
        self.initial_conditions_example = (0.0, 0.0, 0.0, 10.0, 0.0, 0.0)  # I1, I2, VI, Wps, Wk, Wrc
        self.physical_parameters_example = (
            2.0, 2.0,       # L1, L2 (H)
            0.02, 0.02,     # RA1, RA2 (Ohms)
            0.005, 0.005,   # RI1, RI2 (Ohms)
            1000.0,         # CI (F)
            0.1,            # SigmaI (S)
            10.0,           # tauE (s)
            20.0,           # tauk (s)
            60.0            # taurc (s)
        )
        self.t_span_example = (0, 200)
        self.t_eval_example = np.linspace(self.t_span_example[0], self.t_span_example[1], 100)

    def test_runs_with_example_parameters(self):
        """Test if the model runs with the example parameters from windmi_model.py."""
        try:
            solution = solve_windmi_model(
                self.initial_conditions_example,
                self.physical_parameters_example,
                self.t_span_example,
                example_p_in,
                example_p_out_original, # Use the original example_p_out
                t_eval=self.t_eval_example
            )
            self.assertIsNotNone(solution)
            self.assertTrue(solution.success, "Solver should indicate success.")
            self.assertEqual(len(solution.t), len(self.t_eval_example))
        except Exception as e:
            self.fail(f"solve_windmi_model failed with example parameters: {e}")

    def test_zero_input_scenario(self):
        """Test W_ps decay and other variables with zero P_in."""
        initial_conditions_zero_input = (0.0, 0.0, 0.0, 10.0, 0.0, 0.0) # Wps_0 = 10, others 0

        def p_in_zero(t):
            return 0.0

        def p_out_simple_dissipation(t, y_state_vars):
            Wps = y_state_vars[3]
            return 0.1 * Wps # P_out = 0.1 * W_ps

        solution = solve_windmi_model(
            initial_conditions_zero_input,
            self.physical_parameters_example,
            self.t_span_example,
            p_in_zero,
            p_out_simple_dissipation,
            t_eval=self.t_eval_example,
            method='LSODA'  # Specify solver method
        )

        self.assertTrue(solution.success)
        Wps_solution = solution.y[3]

        # Check if W_ps is decaying
        # Allow for small numerical fluctuations by using a small tolerance (atol)
        # Increased tolerance further, as values might be very close to zero at the end
        self.assertTrue(np.all(np.diff(Wps_solution) <= 1e-6), "W_ps should be non-increasing (within tolerance).")
        self.assertLess(Wps_solution[-1], Wps_solution[0], "W_ps should decay from its initial value.")

        # Check other variables (I1, I2, VI) remain small
        # Given P_in is zero and P_out only depends on Wps, I1, I2, VI should not be significantly driven.
        # They might have some transient behavior due to initial conditions of V_I or coupling,
        # but should not grow indefinitely.
        I1_solution = solution.y[0]
        I2_solution = solution.y[1]
        VI_solution = solution.y[2]

        # Allow for some small numerical noise or minor transient effects
        self.assertTrue(np.allclose(I1_solution, 0.0, atol=1e-3), "I1 should remain close to zero.")
        self.assertTrue(np.allclose(I2_solution, 0.0, atol=1e-3), "I2 should remain close to zero.")
        self.assertTrue(np.allclose(VI_solution, 0.0, atol=1e-3), "VI should remain close to zero.")

        # Wk and Wrc are driven by P_out, which depends on Wps. So they will not be zero.
        # Wk = P_out - Wk/tau_k. Wrc = Wk/tau_k - Wrc/tau_rc
        # As Wps decays, P_out decays, and thus Wk and Wrc should also eventually decay.
        Wk_solution = solution.y[4]
        Wrc_solution = solution.y[5]
        # Check they are not negative
        self.assertTrue(np.all(Wk_solution >= -1e-6), "Wk should be non-negative.") # Allow for small numerical errors
        self.assertTrue(np.all(Wrc_solution >= -1e-6), "Wrc should be non-negative.")


    def test_output_shape_and_type(self):
        """Test the shape and type of the output from solve_windmi_model."""
        solution = solve_windmi_model(
            self.initial_conditions_example,
            self.physical_parameters_example,
            self.t_span_example,
            example_p_in,
            example_p_out_original,
            t_eval=self.t_eval_example
        )
        self.assertTrue(solution.success)
        self.assertIsInstance(solution.t, np.ndarray, "sol.t should be a NumPy array.")
        self.assertIsInstance(solution.y, np.ndarray, "sol.y should be a NumPy array.")
        self.assertEqual(solution.y.shape[0], 6, "sol.y should have 6 rows (for 6 state variables).")
        self.assertEqual(solution.y.shape[1], len(self.t_eval_example), "sol.y columns should match t_eval length.")
        self.assertEqual(solution.t.shape[0], len(self.t_eval_example), "sol.t length should match t_eval length.")

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)

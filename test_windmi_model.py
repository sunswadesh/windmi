"""
Unit tests for the WINDMI model implementation in `windmi_model.py`.

This test suite uses the `unittest` framework to verify the correctness
of the `solve_windmi_model` function. It includes tests for general
behavior, specific scenarios like zero input, expected output shape/type,
and validation of the physics-based P_in and P_out calculation formulas.
"""
import unittest
import numpy as np
from windmi_model import solve_windmi_model, heaviside # Import heaviside

class TestWindmiModel(unittest.TestCase):
    """
    Test class for the WINDMI model solver, focusing on the version
    with physics-based P_in and P_out calculations.

    Contains various test methods to ensure the reliability of the
    `solve_windmi_model` function under different conditions.
    """

    def setUp(self):
        """
        Set up common parameters and example inputs for test methods.

        This includes initial conditions, core physical parameters,
        example solar wind parameters (`solar_wind_params_example`) for P_in,
        and example dissipation parameters (`p_out_params_example`) for P_out.
        """
        # Initial conditions: (I1, I2, VI, Wps, Wk, Wrc)
        self.initial_conditions_example = (0.1, 0.1, 1.0, 10.0, 0.1, 0.1)

        # Core physical parameters: (L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc)
        self.physical_parameters_example = (
            2.0, 2.0, 0.02, 0.02, 0.005, 0.005, 1000.0, 0.1, 10.0, 20.0, 60.0
        )

        # Solar wind parameters for P_in
        self.v_sw_fixed = 400.0  # km/s
        self.Bs_imf_fixed = 5.0    # nT (southward component)
        def example_v_sw_func(t): return self.v_sw_fixed
        def example_Bs_imf_func(t): return self.Bs_imf_fixed

        self.solar_wind_params_example = {
            'v_sw_func': example_v_sw_func,
            'Bs_imf_func': example_Bs_imf_func,
            'C_A': 1e-7,  # Placeholder
            'A_eff': (10 * 6371e3)**2  # Placeholder (10*R_E)^2
        }

        # P_out parameters
        self.p_out_params_example = {
            'I_c': 1.0,   # Amperes
            'tau_0': 5.0, # seconds
            'alpha': 2.0  # dimensionless
        }

        self.t_span_example = (0, 50) # Shorter time span for faster tests
        self.t_eval_example = np.linspace(self.t_span_example[0], self.t_span_example[1], 50)


    def test_runs_with_example_parameters(self):
        """
        Test if `solve_windmi_model` runs successfully with the example
        parameters defined in `setUp` (new P_in/P_out physics).
        """
        try:
            solution = solve_windmi_model(
                self.initial_conditions_example,
                self.physical_parameters_example,
                self.solar_wind_params_example,
                self.p_out_params_example,
                self.t_span_example,
                t_eval=self.t_eval_example,
                method='LSODA'
            )
            self.assertIsNotNone(solution)
            self.assertTrue(solution.success, f"Solver failed: {solution.message}")
            self.assertEqual(len(solution.t), len(self.t_eval_example))
        except Exception as e:
            self.fail(f"solve_windmi_model failed with example parameters: {e}")

    def test_zero_input_scenario(self):
        """
        Test W_ps decay when P_in is forced to zero (by setting C_A=0)
        and P_out is made negligible (by setting a very high I_c).
        W_ps should decay primarily due to the -W_ps/tau_E term.
        """
        initial_conditions_zero_input = (0.01, 0.01, 0.1, 10.0, 0.0, 0.0) # Low I1, I2

        # P_in = 0 by setting C_A = 0
        sw_params_zero_p_in = self.solar_wind_params_example.copy()
        sw_params_zero_p_in['C_A'] = 0.0

        # P_out effectively zero by setting I_c very high
        p_out_params_zero_p_out = self.p_out_params_example.copy()
        p_out_params_zero_p_out['I_c'] = 1000.0 # Much larger than initial I1

        solution = solve_windmi_model(
            initial_conditions_zero_input,
            self.physical_parameters_example,
            sw_params_zero_p_in,
            p_out_params_zero_p_out,
            self.t_span_example,
            t_eval=self.t_eval_example,
            method='LSODA'
        )

        self.assertTrue(solution.success, f"Solver failed: {solution.message}")
        Wps_solution = solution.y[3]
        tauE = self.physical_parameters_example[8]

        # With P_in=0 and P_out=0 (due to high I_c), dWps/dt = -Wps/tauE
        # So Wps should decay exponentially.
        self.assertTrue(np.all(np.diff(Wps_solution) < 1e-6), # Allow small numerical errors for "non-increasing"
                        "W_ps should be generally non-increasing.")
        self.assertLess(Wps_solution[-1], Wps_solution[0],
                        "W_ps should decay from its initial value.")

        # Check if P_out was indeed zero or very close to it
        self.assertTrue(hasattr(solution, 'p_out_values'), "Solution object should have p_out_values")
        self.assertTrue(np.allclose(solution.p_out_values, 0.0, atol=1e-9), "P_out should be effectively zero.")


    def test_output_shape_and_type(self):
        """
        Test the shape and type of the output from `solve_windmi_model`
        using the new P_in/P_out physics.
        Ensures sol.y, sol.t, sol.p_in_values, sol.p_out_values are correct.
        """
        solution = solve_windmi_model(
            self.initial_conditions_example,
            self.physical_parameters_example,
            self.solar_wind_params_example,
            self.p_out_params_example,
            self.t_span_example,
            t_eval=self.t_eval_example,
            method='LSODA'
        )
        self.assertTrue(solution.success, f"Solver failed: {solution.message}")
        self.assertIsInstance(solution.t, np.ndarray, "sol.t should be a NumPy array.")
        self.assertIsInstance(solution.y, np.ndarray, "sol.y should be a NumPy array.")
        self.assertEqual(solution.y.shape[0], 6, "sol.y should have 6 rows.")
        self.assertEqual(solution.y.shape[1], len(self.t_eval_example), "sol.y columns should match t_eval.")
        self.assertEqual(solution.t.shape[0], len(self.t_eval_example), "sol.t length should match t_eval.")
        self.assertTrue(hasattr(solution, 'p_in_values'), "Solution object should have p_in_values")
        self.assertTrue(hasattr(solution, 'p_out_values'), "Solution object should have p_out_values")
        self.assertEqual(len(solution.p_in_values), len(self.t_eval_example), "p_in_values length mismatch")
        self.assertEqual(len(solution.p_out_values), len(self.t_eval_example), "p_out_values length mismatch")


    def test_p_in_calculation(self):
        """
        Test if P_in is calculated correctly according to the new formula
        P_in(t) = C_A * v_sw(t) * Bs_imf(t) * A_eff, using constant inputs.
        """
        # Use fixed solar wind parameters from setUp
        v_sw = self.solar_wind_params_example['v_sw_func'](0)
        Bs_imf = self.solar_wind_params_example['Bs_imf_func'](0)
        C_A = self.solar_wind_params_example['C_A']
        A_eff = self.solar_wind_params_example['A_eff']
        expected_p_in = C_A * v_sw * Bs_imf * A_eff

        solution = solve_windmi_model(
            self.initial_conditions_example,
            self.physical_parameters_example,
            self.solar_wind_params_example, # Uses the functions returning fixed values
            self.p_out_params_example,
            self.t_span_example,
            t_eval=self.t_eval_example,
            method='LSODA'
        )
        self.assertTrue(solution.success, f"Solver failed: {solution.message}")
        self.assertTrue(hasattr(solution, 'p_in_values'), "Solution object should have p_in_values")
        self.assertTrue(np.allclose(solution.p_in_values, expected_p_in, rtol=1e-5),
                        f"P_in values {solution.p_in_values} not close to expected {expected_p_in}")

    def test_p_out_calculation_below_critical_current(self):
        """
        Test if P_out (new formula) is zero when magnetotail current I1
        is below the critical current I_c.
        """
        ic_test = (0.01, 0.01, 0.1, 10.0, 0.1, 0.1) # Low I1

        p_out_params_high_Ic = self.p_out_params_example.copy()
        p_out_params_high_Ic['I_c'] = 100.0 # I_c is much higher than I1

        solution = solve_windmi_model(
            ic_test,
            self.physical_parameters_example,
            self.solar_wind_params_example,
            p_out_params_high_Ic,
            self.t_span_example,
            t_eval=self.t_eval_example,
            method='LSODA'
        )
        self.assertTrue(solution.success, f"Solver failed: {solution.message}")
        self.assertTrue(hasattr(solution, 'p_out_values'), "Solution object should have p_out_values")
        self.assertTrue(np.all(solution.p_out_values == 0.0),
                        f"P_out values {solution.p_out_values} should be all zero when I1 < I_c.")

    def test_p_out_calculation_above_critical_current(self):
        """
        Test P_out calculation using the new formula when magnetotail
        current I1 is above the critical current I_c.
        Verifies calculated P_out against expected values for initial time steps.
        """
        # Ensure I1 starts above I_c and stays there for a bit, or is driven there
        I1_initial = 2.0
        Ic_val = 1.0
        Wps_initial = 20.0
        tau0_val = self.p_out_params_example['tau_0']
        alpha_val = self.p_out_params_example['alpha']

        ic_test = (I1_initial, 0.1, 1.0, Wps_initial, 0.1, 0.1)

        p_out_params_active = self.p_out_params_example.copy()
        p_out_params_active['I_c'] = Ic_val
        p_out_params_active['alpha'] = alpha_val
        p_out_params_active['tau_0'] = tau0_val

        # For this test, make P_in small or zero so it doesn't overwhelm Wps changes quickly
        sw_params_low_pin = self.solar_wind_params_example.copy()
        sw_params_low_pin['C_A'] = 1e-10 # Very small P_in

        # Short time span to check initial P_out
        t_span_short = (0, 2.0)
        t_eval_short = np.linspace(t_span_short[0], t_span_short[1], 10)

        solution = solve_windmi_model(
            ic_test,
            self.physical_parameters_example,
            sw_params_low_pin,
            p_out_params_active,
            t_span_short,
            t_eval=t_eval_short,
            method='LSODA'
        )
        self.assertTrue(solution.success, f"Solver failed: {solution.message}")
        self.assertTrue(hasattr(solution, 'p_out_values'), "Solution object should have p_out_values")

        # Check P_out at the first few points
        for i in range(min(3, len(solution.t))): # Check first few points
            I1_sol = solution.y[0, i]
            Wps_sol = solution.y[3, i]

            expected_p_out_val = 0.0
            if Ic_val > 0:
                expected_p_out_val = (Wps_sol / tau0_val) * \
                                   heaviside(I1_sol - Ic_val) * \
                                   ((I1_sol - Ic_val) / Ic_val)**alpha_val

            self.assertAlmostEqual(solution.p_out_values[i], expected_p_out_val, places=5,
                                   msg=f"P_out at t={solution.t[i]} mismatch. Got {solution.p_out_values[i]}, expected {expected_p_out_val}")
            if I1_initial > Ic_val: # Expect P_out to be non-zero if Wps_initial > 0
                 self.assertGreater(solution.p_out_values[0], 0, "P_out should be > 0 initially if I1 > Ic and Wps > 0")


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)

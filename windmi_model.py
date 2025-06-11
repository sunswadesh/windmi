"""
Implements the Wind-Magnetosphere-Ionosphere (WINDMI) model.

The WINDMI model is a physics-based nonlinear dynamical model for the
nightside magnetosphere, including the magnetotail and the ring current.
It describes the global dynamics of the magnetosphere-ionosphere system
driven by solar wind input. The model uses an equivalent solar wind
voltage source as input to a system of six ordinary differential equations (ODEs)
representing the state of the magnetotail lobes, central plasma sheet,
and ring current.

This implementation is based on the work by:
- Doxas, I., J. R. Spencer, J. W. Freeman, and P. H. Reiff (2004),
  A new nonlinear dynamical model of the magnetotail, J. Geophys. Res.,
  109, A10209, doi:10.1029/2003JA010225.
- Spencer, J. R., I. Doxas, J. W. Freeman, and P. H. Reiff (2007),
  Nonlinear dynamics of the magnetotail: The WINDMI model,
  in Modeling of the Space Weather Environment, edited by J. M. Schmidtke
  and D. L. Gallagher, pp. 199–224, AGU, Washington, D. C.
- Patra, A., I. Doxas, and D. N. Baker (2011), Data assimilation in the
  WINDMI model of the magnetosphere, Space Weather, 9, S04006,
  doi:10.1029/2010SW000618.

The model equations for the state variables (I_1, I_2, V_I, W_ps, W_k, W_rc) are:
1. dI_1/dt = (V_I - (R_A1 + R_I1) * I_1) / L_1
2. dI_2/dt = (V_I - (R_A2 + R_I2) * I_2) / L_2
3. dV_I/dt = (I_1 + I_2 - Σ_I * V_I) / C_I
4. dW_ps/dt = P_in - W_ps/τ_E - P_out
5. dW_k/dt = P_out - W_k/τ_k
6. dW_rc/dt = W_k/τ_k - W_rc/τ_rc

The power input P_in is calculated as:
P_in(t) = C_A * v_sw(t) * Bs_imf(t) * A_eff
where v_sw(t) is the solar wind speed, Bs_imf(t) is the magnitude of the
southward component of the Interplanetary Magnetic Field (IMF), C_A is a
coupling coefficient, and A_eff is an effective area.

The power output P_out (energy dissipation) is calculated as:
P_out(t, y) = (W_ps / τ_0) * H(I_1 - I_c) * ((I_1 - I_c) / I_c)^α
where W_ps is the plasma sheet energy, I_1 is the dawn sector current,
H is the Heaviside step function, and I_c, τ_0, α are parameters
controlling the dissipation process.

Core physical parameters (L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc)
and parameters for P_in/P_out calculations are passed to the solver function.
"""
import numpy as np
from scipy.integrate import solve_ivp

def heaviside(x):
    """
    Heaviside step function.

    Parameters
    ----------
    x : float or int
        Input value.

    Returns
    -------
    int
        1 if x >= 0, else 0.
    """
    return 1 if x >= 0 else 0

def solve_windmi_model(initial_conditions, physical_parameters,
                       solar_wind_params, p_out_params,
                       t_span, t_eval=None, **kwargs):
    """
    Solves the WINDMI model Ordinary Differential Equation (ODE) system
    with updated P_in and P_out formulations.

    Parameters
    ----------
    initial_conditions : tuple or array_like
        Initial values for the six state variables, in the order:
        - I_1: Current in the dawn sector of the magnetotail (Amperes)
        - I_2: Current in the dusk sector of the magnetotail (Amperes)
        - V_I: Cross-polar cap potential (Volts)
        - W_ps: Energy stored in the pressure gradient of the plasma sheet (Joules)
        - W_k: Bulk kinetic energy of plasma flow in the plasma sheet (Joules)
        - W_rc: Energy stored in the ring current (Joules)
    physical_parameters : tuple or array_like
        Values for the eleven core physical parameters of the model, in the order:
        - L_1: Inductance of the dawn sector magnetotail (Henrys)
        - L_2: Inductance of the dusk sector magnetotail (Henrys)
        - R_A1: Resistance of the dawn auroral ionosphere (Ohms)
        - R_A2: Resistance of the dusk auroral ionosphere (Ohms)
        - R_I1: Resistance of the dawn sector magnetotail current path (Ohms)
        - R_I2: Resistance of the dusk sector magnetotail current path (Ohms)
        - C_I: Capacitance of the ionosphere (Farads)
        - Σ_I: Conductance of the ionosphere (Siemens)
        - τ_E: Energy unloading timescale for plasma sheet pressure (seconds)
        - τ_k: Timescale for kinetic energy transfer to ring current or dissipation (seconds)
        - τ_rc: Decay timescale for ring current energy (seconds)
    solar_wind_params : tuple or dict
        Parameters required for the P_in calculation.
        If a tuple, expected order is `(v_sw_func, Bs_imf_func, C_A, A_eff)`.
        If a dictionary, expected keys are `'v_sw_func'`, `'Bs_imf_func'`, `'C_A'`, `'A_eff'`.
            - `v_sw_func` : callable
                Function `v_sw(t)` that returns solar wind speed at time `t`.
                Example units: km/s.
            - `Bs_imf_func` : callable
                Function `Bs_imf(t)` that returns the magnitude of the southward
                component of the IMF (B_s) at time `t`. This value is typically
                non-negative (B_s = |Bz| if Bz < 0, else 0, or total B_s magnitude).
                Example units: nT.
            - `C_A` : float
                Coupling coefficient for P_in. Its units must be consistent with
                `v_sw`, `Bs_imf`, and `A_eff` to result in P_in in Watts.
            - `A_eff` : float
                Effective area for solar wind energy input. Example units: m^2.
    p_out_params : tuple or dict
        Parameters required for the P_out calculation.
        If a tuple, expected order is `(I_c, tau_0, alpha)`.
        If a dictionary, expected keys are `'I_c'`, `'tau_0'`, `'alpha'`.
            - `I_c` : float
                Critical current for the activation of P_out. Units: Amperes.
            - `tau_0` : float
                Characteristic timescale for energy dissipation in P_out. Units: seconds.
            - `alpha` : float
                Exponent in the P_out formula, controlling the non-linearity of
                the dissipation. Dimensionless.
    t_span : tuple
        Time interval for the simulation `(t_start, t_end)`. Units: seconds.
    t_eval : array_like, optional
        Times at which to store the computed solution. Must be sorted and lie
        within `t_span`. If None (default), solution is stored at points chosen
        by the solver.
    **kwargs : dict, optional
        Additional keyword arguments to pass to `scipy.integrate.solve_ivp`
        (e.g., `method`, `rtol`, `atol`).

    Returns
    -------
    scipy.integrate.OdeResult
        An object containing the solution of the ODE system. Key attributes include:
        - `t` : ndarray, shape (n_points,)
            Time points.
        - `y` : ndarray, shape (6, n_points)
            Values of the state variables at the time points. Each row
            corresponds to a state variable in the order defined in
            `initial_conditions` (e.g., `y[0]` is I_1(t), `y[1]` is I_2(t), etc.).
        - `sol` : OdeSolution (if `dense_output=True` was used, which is default here)
            A callable function to evaluate the solution at any time `t` within `t_span`.
        - `success` : bool
            True if the solver reached the end of the integration interval.
        - `p_in_values` : ndarray, shape (n_points,)
            Calculated P_in values at each time point in `sol.t`. Only if `t_eval` is provided.
        - `p_out_values` : ndarray, shape (n_points,)
            Calculated P_out values at each time point in `sol.t`. Only if `t_eval` is provided.
    """
    I1_0, I2_0, VI_0, Wps_0, Wk_0, Wrc_0 = initial_conditions
    L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc = physical_parameters

    # Unpack solar_wind_params and p_out_params for use in windmi_odes and P_in/P_out calculation
    if isinstance(solar_wind_params, dict):
        v_sw_func_local = solar_wind_params['v_sw_func']
        Bs_imf_func_local = solar_wind_params['Bs_imf_func']
        C_A_local = solar_wind_params['C_A']
        A_eff_local = solar_wind_params['A_eff']
    else:
        v_sw_func_local, Bs_imf_func_local, C_A_local, A_eff_local = solar_wind_params

    if isinstance(p_out_params, dict):
        I_c_local = p_out_params['I_c']
        tau_0_local = p_out_params['tau_0']
        alpha_local = p_out_params['alpha']
    else:
        I_c_local, tau_0_local, alpha_local = p_out_params


    def windmi_odes(t, y):
        """
        Defines the system of ODEs for the WINDMI model.
        Uses parameters from the outer scope (L1..taurc, v_sw_func_local.., I_c_local..).
        """
        I1, I2, VI, Wps, Wk, Wrc = y

        # Calculate P_in using the new formula
        val_v_sw = v_sw_func_local(t)
        val_Bs_imf = Bs_imf_func_local(t) # Assumed to be southward component magnitude
        current_p_in = C_A_local * val_v_sw * val_Bs_imf * A_eff_local

        # Calculate P_out using the new formula
        if I_c_local > 0: # Ensure I_c is positive to avoid division by zero or math errors
            # Power term is active only if I1 > I_c
            # The term ((I1 - I_c) / I_c)**alpha can be complex if (I1 - I_c) < 0 and alpha is not an integer.
            # Heaviside function ensures this term is zero if I1 - I_c < 0.
            # If I1 - I_c == 0, P_out is 0.
            # If I1 - I_c > 0, then (I1 - I_c)/I_c is positive.
            if I1 > I_c_local:
                 trigger_term = ((I1 - I_c_local) / I_c_local)**alpha_local
            else:
                 trigger_term = 0.0 # Handles I1 <= I_c due to Heaviside logic
            current_p_out = (Wps / tau_0_local) * trigger_term * heaviside(I1 - I_c_local)

        else: # If I_c is not positive, P_out is zero.
            current_p_out = 0.0

        # The WINDMI ODEs:
        dI1_dt = (VI - (RA1 + RI1) * I1) / L1  # Eq. 1
        dI2_dt = (VI - (RA2 + RI2) * I2) / L2  # Eq. 2
        dVI_dt = (I1 + I2 - SigmaI * VI) / CI   # Eq. 3
        dWps_dt = current_p_in - Wps / tauE - current_p_out # Eq. 4
        dWk_dt = current_p_out - Wk / tauk      # Eq. 5
        dWrc_dt = Wk / tauk - Wrc / taurc     # Eq. 6

        return [dI1_dt, dI2_dt, dVI_dt, dWps_dt, dWk_dt, dWrc_dt]

    sol = solve_ivp(
        windmi_odes,
        t_span,
        [I1_0, I2_0, VI_0, Wps_0, Wk_0, Wrc_0], # Initial state vector
        t_eval=t_eval,
        dense_output=True,  # Enable dense output for smoother plots if needed
        **kwargs
    )

    # Store calculated P_in and P_out for analysis if needed
    # This requires re-evaluating them at solution time points.
    # Note: This is an approximation if the solver used adaptive time steps
    # different from t_eval. For precise values at each internal step,
    # modification of the ODE solver loop would be needed, or by using sol.sol(t).
    # The dense_output=True (default in solve_ivp used here) ensures sol.sol is available.
    if sol.success and t_eval is not None:
        sol.p_in_values = np.array([
            C_A_local * v_sw_func_local(t) * Bs_imf_func_local(t) * A_eff_local for t in sol.t
        ])
        sol.p_out_values = np.zeros_like(sol.t)
        for i, t_val in enumerate(sol.t):
            y_at_t = sol.sol(t_val) # Get state variables at time t_val
            I1_at_t = y_at_t[0]
            Wps_at_t = y_at_t[3]
            if I_c_local > 0:
                if I1_at_t > I_c_local: # Calculate P_out only if I1 > I_c
                    trigger = ((I1_at_t - I_c_local) / I_c_local)**alpha_local
                    sol.p_out_values[i] = (Wps_at_t / tau_0_local) * trigger * heaviside(I1_at_t - I_c_local) # heaviside for safety
                else:
                    sol.p_out_values[i] = 0.0
            else: # If I_c is not positive, P_out is zero.
                sol.p_out_values[i] = 0.0
    else: # If solver failed or t_eval was None, return empty arrays for p_in/p_out
        sol.p_in_values = np.array([])
        sol.p_out_values = np.array([])

    return sol


if __name__ == '__main__':
    # Example Usage for the updated WINDMI model with physics-based P_in and P_out

    # 1. Define Initial Conditions for the 6 state variables:
    # (I1, I2, VI, Wps, Wk, Wrc)
    # Units: Amperes, Amperes, Volts, Joules, Joules, Joules
    ic = (0.1, 0.1, 1.0, 10.0, 0.1, 0.1)  # Example values

    # 2. Define Core Physical Parameters for the system (11 parameters):
    # (L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc)
    # Units: H, H, Ohm, Ohm, Ohm, Ohm, F, S, s, s, s respectively
    physical_params = (
        2.0, 2.0,       # L1, L2 (Henrys)
        0.02, 0.02,     # RA1, RA2 (Ohms) - Auroral ionosphere resistances
        0.005, 0.005,   # RI1, RI2 (Ohms) - Magnetotail current path resistances
        1000.0,         # CI (Farads) - Ionospheric capacitance
        0.1,            # SigmaI (Siemens) - Ionospheric conductance
        10.0,           # tauE (seconds) - Plasma sheet energy unloading timescale
        20.0,           # tauk (seconds) - Kinetic energy transfer timescale
        60.0            # taurc (seconds) - Ring current decay timescale
    )

    # 3. Define Solar Wind Parameters for P_in calculation
    # These include functions for time-varying solar wind speed and IMF Bs,
    # a coupling coefficient C_A, and an effective area A_eff.

    # Example: Define a function for solar wind speed v_sw(t)
    def example_v_sw(t): # Example returns speed in km/s
        return 400.0 # Constant solar wind speed of 400 km/s

    # Example: Define a function for IMF Southward Component Bs_imf(t)
    def example_Bs_imf(t): # Example returns Bs in nT (should be positive for southward)
        # Simulate a step increase in Bs after 50 seconds
        if t < 50:
            return 1.0 # nT
        else:
            return 5.0 # nT

    # Define coupling coefficient C_A and effective area A_eff.
    # Note: The physical units of C_A must be chosen carefully to ensure P_in is in Watts,
    #       depending on the units of v_sw (e.g., m/s vs km/s) and Bs_imf (e.g., T vs nT).
    #       P_in = C_A * v_sw * Bs_imf * A_eff.
    #       If v_sw is in m/s, Bs_imf in Tesla, A_eff in m^2, then C_A is dimensionless
    #       (or related to Poynting flux conversion).
    #       Here, v_sw(km/s) * 1e3 -> m/s. Bs_imf(nT) * 1e-9 -> T.
    #       So, P_in = C_A_val * (v_sw * 1e3) * (Bs_imf * 1e-9) * A_eff_val
    #       To make C_A_val closer to 1, one might absorb 1e3*1e-9 into it if inputs are always km/s and nT.
    #       For this example, we keep C_A as a direct multiplier.
    C_A_val = 1.0 # Example: A simplified coefficient. For physical realism, this needs calibration.
                  # (e.g. related to 1/mu_0 if P_in ~ Poynting flux * A_eff)
    A_eff_val = (10 * 6371e3)**2 # Effective area example: (10 * Earth_Radius)^2 in m^2

    # Store solar wind parameters in a dictionary
    sw_params = {
        'v_sw_func': example_v_sw,      # Callable function for solar wind speed
        'Bs_imf_func': example_Bs_imf,  # Callable function for IMF Bs
        'C_A': C_A_val,                 # Coupling coefficient
        'A_eff': A_eff_val              # Effective area
    }

    # 4. Define Parameters for P_out calculation
    # These include critical current I_c, dissipation timescale tau_0, and exponent alpha.
    p_out_params_vals = {
        'I_c': 1.0,     # Critical current for P_out activation (Amperes)
        'tau_0': 5.0,   # Characteristic dissipation timescale (seconds)
        'alpha': 2.0    # Exponent in P_out formula (dimensionless)
    }

    # 5. Define Time Span for Simulation (in seconds)
    t_start = 0
    t_end = 200  # seconds
    tspan = (t_start, t_end)
    t_points = np.linspace(t_start, t_end, 500) # evaluation points

    # 6. Run Simulation
    print("Running WINDMI model simulation with updated P_in/P_out formulas...")
    solution = solve_windmi_model(
        ic, physical_params, sw_params, p_out_params_vals,
        tspan, t_eval=t_points, method='LSODA' # Using LSODA for stiffness
    )
    print("Simulation complete.")

    # 7. Accessing and Plotting Results
    if solution.success:
        time = solution.t
        I1_sol = solution.y[0]
        I2_sol = solution.y[1]
        VI_sol = solution.y[2]
        Wps_sol = solution.y[3]
        Wk_sol = solution.y[4]
        Wrc_sol = solution.y[5]

        p_in_calculated = solution.p_in_values
        p_out_calculated = solution.p_out_values

        print(f"Number of time points: {len(time)}")
        print(f"Final Wps: {Wps_sol[-1] if len(Wps_sol) > 0 else 'N/A'}")

        try:
            import matplotlib.pyplot as plt
            print("Plotting results...")
            plt.figure(figsize=(12, 12)) # Adjusted figure size

            plt.subplot(4, 2, 1)
            plt.plot(time, I1_sol, label='I_1(t)')
            plt.plot(time, I2_sol, label='I_2(t)', linestyle='--')
            plt.xlabel('Time (s)'); plt.ylabel('Currents (A)'); plt.legend(); plt.title('Magnetotail Currents')

            plt.subplot(4, 2, 2)
            plt.plot(time, VI_sol, label='V_I(t)')
            plt.xlabel('Time (s)'); plt.ylabel('Potential (V)'); plt.legend(); plt.title('Cross-Polar Cap Potential')

            plt.subplot(4, 2, 3)
            plt.plot(time, Wps_sol, label='W_ps(t)')
            plt.xlabel('Time (s)'); plt.ylabel('Energy (J)'); plt.legend(); plt.title('Pressure-Gradient Energy')

            plt.subplot(4, 2, 4)
            plt.plot(time, Wk_sol, label='W_k(t)')
            plt.xlabel('Time (s)'); plt.ylabel('Energy (J)'); plt.legend(); plt.title('Bulk Kinetic Energy')

            plt.subplot(4, 2, 5)
            plt.plot(time, Wrc_sol, label='W_rc(t)')
            plt.xlabel('Time (s)'); plt.ylabel('Energy (J)'); plt.legend(); plt.title('Ring Current Energy')

            plt.subplot(4, 2, 6)
            if p_in_calculated.size > 0:
                 plt.plot(time, p_in_calculated, label='P_in(t) Calculated', color='green')
            if p_out_calculated.size > 0:
                 plt.plot(time, p_out_calculated, label='P_out(t) Calculated', color='red', linestyle='--')
            plt.xlabel('Time (s)'); plt.ylabel('Power (W)'); plt.legend(); plt.title('Input and Output Power')

            plt.subplot(4, 2, 7)
            plt.plot(time, [sw_params['v_sw_func'](t) for t in time], label='v_sw(t)')
            plt.xlabel('Time (s)'); plt.ylabel('Solar Wind Speed (km/s)'); plt.legend(); plt.title('Solar Wind Speed Input')

            plt.subplot(4, 2, 8)
            plt.plot(time, [sw_params['Bs_imf_func'](t) for t in time], label='Bs_imf(t)')
            plt.xlabel('Time (s)'); plt.ylabel('IMF Bs (nT)'); plt.legend(); plt.title('IMF Bs Input')


            plt.tight_layout()
            plt.show()
            print("Plotting complete. If plots are not showing, ensure you have a GUI environment.")
        except ImportError:
            print("Matplotlib not found. Skipping plotting.")
        except Exception as e:
            print(f"An error occurred during plotting: {e}")
    else:
        print(f"Solver failed: {solution.message}")

    print("\nNote: Example parameters for P_in (C_A, A_eff) and solar wind functions are illustrative.")
    print("These require careful calibration for physically meaningful results.")

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

The model equations are:
1. dI_1/dt = (V_I - (R_A1 + R_I1) * I_1) / L_1
2. dI_2/dt = (V_I - (R_A2 + R_I2) * I_2) / L_2
3. dV_I/dt = (I_1 + I_2 - Σ_I * V_I) / C_I
4. dW_ps/dt = P_in - W_ps/τ_E - P_out
5. dW_k/dt = P_out - W_k/τ_k
6. dW_rc/dt = W_k/τ_k - W_rc/τ_rc
"""
import numpy as np
from scipy.integrate import solve_ivp

def solve_windmi_model(initial_conditions, physical_parameters, t_span, p_in, p_out, t_eval=None, **kwargs):
    """
    Solves the WINDMI model Ordinary Differential Equation (ODE) system.

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
        Values for the eleven physical parameters of the model, in the order:
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
    t_span : tuple
        Time interval for the simulation (t_start, t_end) in seconds.
    p_in : float or callable
        Solar wind input power (Watts). This is often calculated from solar wind
        parameters like velocity (v) and the southward component of the
        Interplanetary Magnetic Field (IMF B_s), typically using a rectified
        vBs formula (e.g., P_in is proportional to v * B_s^2 * sin^4(theta/2),
        where theta is the IMF clock angle).
        If callable, it should be a function of time `t`, `p_in(t)`.
    p_out : float or callable
        Power output term (Watts), representing energy transfer from the plasma
        sheet, primarily to the ring current and auroral dissipation.
        If callable, it should be a function of time `t` and the state
        variables `y = (I_1, I_2, V_I, W_ps, W_k, W_rc)`, i.e., `p_out(t, y)`.
    t_eval : array_like, optional
        Times at which to store the computed solution. Must be sorted and lie
        within `t_span`. If None (default), the points are chosen by the solver.
    **kwargs : dict, optional
        Additional keyword arguments to pass to `scipy.integrate.solve_ivp`
        (e.g., `method`, `rtol`, `atol`).

    Returns
    -------
    scipy.integrate.OdeResult
        An object containing the solution of the ODE system. Key attributes:
        - t : ndarray, shape (n_points,)
            Time points.
        - y : ndarray, shape (6, n_points)
            Values of the state variables at the time points. Each row
            corresponds to a state variable in the order defined in
            `initial_conditions`:
            - y[0]: I_1(t)
            - y[1]: I_2(t)
            - y[2]: V_I(t)
            - y[3]: W_ps(t)
            - y[4]: W_k(t)
            - y[5]: W_rc(t)
        - sol : OdeSolution (if `dense_output=True`)
            A callable function to evaluate the solution at any time t within t_span.
        - success : bool
            True if the solver reached the end of the integration interval.
    """
    I1_0, I2_0, VI_0, Wps_0, Wk_0, Wrc_0 = initial_conditions
    L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc = physical_parameters

    def windmi_odes(t, y):
        """Defines the system of ODEs for the WINDMI model."""
        I1, I2, VI, Wps, Wk, Wrc = y

        # Evaluate P_in and P_out if they are functions
        current_p_in = p_in(t) if callable(p_in) else p_in
        current_p_out = p_out(t, y) if callable(p_out) else p_out

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

    return sol

def example_p_in(t):
    """
    Example solar wind input power (P_in) function for demonstration.

    This function simulates a gradual increase and then decrease in power input.
    It is purely illustrative and not based on actual solar wind data.
    A realistic P_in function would typically derive power from solar wind
    parameters (e.g., velocity v, IMF B_s).

    Parameters
    ----------
    t : float
        Time in seconds.

    Returns
    -------
    float
        Calculated input power P_in (Watts) at time t.
    """
    if t < 10: # s
        return 0.0 # Watts
    elif 10 <= t < 100: # s
        # Simulate a gradual increase (arbitrary units for illustration)
        return 10.0 * (1 - np.exp(-(t - 10.0) / 20.0)) # Watts
    else: # s
        # Simulate a gradual decrease (arbitrary units for illustration)
        return 10.0 * np.exp(-(t - 100.0) / 50.0) # Watts

def example_p_out(t, y_state_vars):
    """
    Example power output (P_out) function for demonstration.

    This function provides a highly simplified placeholder for P_out, making it
    dependent on the plasma sheet energy (W_ps) and ionospheric potential (V_I).
    A more realistic P_out would be based on physical models of energy
    dissipation, such as auroral precipitation (e.g., related to AL index)
    and Joule heating.

    Parameters
    ----------
    t : float
        Time in seconds.
    y_state_vars : array_like
        Current state variables (I1, I2, VI, Wps, Wk, Wrc).
        Units: Amperes, Amperes, Volts, Joules, Joules, Joules.

    Returns
    -------
    float
        Calculated output power P_out (Watts) at time t for the given state.
    """
    I1, I2, VI, Wps, Wk, Wrc = y_state_vars
    # Example: P_out is proportional to Wps and a saturating function of VI.
    # This is a placeholder and needs to be defined based on physical models.
    k_dissipation = 0.05 # Arbitrary dissipation constant for this example (unitless)
    # Wps (J), VI (V)
    return k_dissipation * Wps * (abs(VI) / (1.0 + abs(VI))) # Watts

if __name__ == '__main__':
    # Example Usage: Demonstrates how to run the model with placeholder P_in and P_out.
    # Define initial conditions (example values)
    # (I1, I2, VI, Wps, Wk, Wrc)
    ic = (0.0, 0.0, 0.0, 10.0, 0.0, 0.0)  # Units: A, A, V, J, J, J

    # Define physical parameters (example values from typical ranges)
    # (L1, L2, RA1, RA2, RI1, RI2, CI, SigmaI, tauE, tauk, taurc)
    pp = (
        2.0, 2.0,       # L1, L2 (Henrys)
        0.02, 0.02,     # RA1, RA2 (Ohms) - Auroral ionosphere resistances
        0.005, 0.005,   # RI1, RI2 (Ohms) - Magnetotail current path resistances
        1000.0,         # CI (Farads) - Ionospheric capacitance
        0.1,            # SigmaI (Siemens) - Ionospheric conductance
        10.0,           # tauE (seconds) - Plasma sheet energy unloading timescale
        20.0,           # tauk (seconds) - Kinetic energy transfer timescale
        60.0            # taurc (seconds) - Ring current decay timescale
    )

    # Define time span for simulation (seconds)
    t_start = 0
    t_end = 200
    tspan = (t_start, t_end)
    t_points = np.linspace(t_start, t_end, 500) # evaluation points

    print("Running WINDMI model simulation with example parameters...")
    solution = solve_windmi_model(ic, pp, tspan, example_p_in, example_p_out, t_eval=t_points)
    print("Simulation complete.")

    # Accessing results:
    # Time points
    time = solution.t
    # State variables (each row is a variable, each column is a time point)
    I1_sol = solution.y[0]
    I2_sol = solution.y[1]
    VI_sol = solution.y[2]
    Wps_sol = solution.y[3]
    Wk_sol = solution.y[4]
    Wrc_sol = solution.y[5]

    print(f"Number of time points: {len(time)}")
    print(f"Final Wps: {Wps_sol[-1]}")

    # Optional: Plotting results (requires matplotlib)
    try:
        import matplotlib.pyplot as plt
        print("Plotting results...")
        plt.figure(figsize=(12, 10))

        plt.subplot(3, 2, 1)
        plt.plot(time, I1_sol, label='I_1(t)')
        plt.plot(time, I2_sol, label='I_2(t)', linestyle='--')
        plt.xlabel('Time (s)')
        plt.ylabel('Currents (A)')
        plt.legend()
        plt.title('Magnetotail Currents')

        plt.subplot(3, 2, 2)
        plt.plot(time, VI_sol, label='V_I(t)')
        plt.xlabel('Time (s)')
        plt.ylabel('Ionospheric Potential (V)')
        plt.legend()
        plt.title('Cross-Polar Cap Potential')

        plt.subplot(3, 2, 3)
        plt.plot(time, Wps_sol, label='W_ps(t)')
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (J)')
        plt.legend()
        plt.title('Pressure-Gradient Energy')

        plt.subplot(3, 2, 4)
        plt.plot(time, Wk_sol, label='W_k(t)')
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (J)')
        plt.legend()
        plt.title('Bulk Kinetic Energy of Plasma Flow')

        plt.subplot(3, 2, 5)
        plt.plot(time, Wrc_sol, label='W_rc(t)')
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (J)')
        plt.legend()
        plt.title('Ring Current Energy')

        # Plot P_in and P_out
        p_in_values = [example_p_in(t_val) for t_val in time]
        p_out_values = [example_p_out(t_val, solution.sol(t_val)) for t_val in time]

        plt.subplot(3, 2, 6)
        plt.plot(time, p_in_values, label='P_in(t)', color='green')
        plt.plot(time, p_out_values, label='P_out(t)', color='red', linestyle='--')
        plt.xlabel('Time (s)')
        plt.ylabel('Power (W)')
        plt.legend()
        plt.title('Input and Output Power')


        plt.tight_layout()
        plt.show()
        print("Plotting complete. If plots are not showing, ensure you have a GUI environment.")
    except ImportError:
        print("Matplotlib not found. Skipping plotting.")
    except Exception as e:
        print(f"An error occurred during plotting: {e}")

    print("\nNote: The P_in and P_out functions in this example are placeholders.")
    print("For a realistic simulation, these need to be defined based on actual solar wind data and physical models of energy dissipation.")

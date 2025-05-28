## WINDMI Model Description

The WINDMI (WINd-Magnetosphere-Ionosphere) model is a scientific tool designed to simulate and understand the complex interactions between the solar wind and Earth's magnetosphere and ionosphere.

**Purpose of the Model:**
The primary purpose of the WINDMI model is to simulate and predict the transfer of energy from the solar wind, through the Earth's magnetosphere, and into the ionosphere. It aims to estimate and forecast space weather phenomena, particularly geomagnetic activity and disturbances.

**Physical System Modeled:**
WINDMI models the coupled solar wind-magnetosphere-ionosphere system. This includes:
*   The **solar wind** as the external driver.
*   The **magnetosphere**, focusing on aspects like the ring current and tail current.
*   The **ionosphere**, including its conductivity and currents (Region I and II).

**Modeling Approach:**
WINDMI is a low-dimensional (low-order) model. It employs an **electric circuitry analogy** to represent the magnetosphere-ionosphere system. This means components like capacitances, resistances, and inductances are used to conceptualize and mathematically describe the system's behavior. The model is governed by a set of **nonlinear ordinary differential equations (ODEs)**. The electric driving voltage applied by the solar wind can be described using various coupling functions, such as the Rectified (Reiff and Luhmann, 1986), Siscoe (Siscoe et al. 2002), or Newell (2007) functions. While the core model is physics-based, modular machine learning tools can be used for parameter optimization.

**Typical Inputs to the Model:**
The model is typically driven by:
*   **Solar wind measurements:** These can be real-time (ACE Real-Time) or archived (ACE Level2) data. Key parameters from the solar wind (like velocity, density, and interplanetary magnetic field components) are used to calculate the driving voltage.

**Typical Outputs or Products of the Model:**
The main outputs and products of the WINDMI model include:
*   **Energy of the ring current:** This is often represented by or used to predict geomagnetic indices like Dst or SymH.
*   **Auroral indices:** Such as AU and AL, which quantify the intensity of auroral electrojets.
*   **Other magnetospheric and ionospheric parameters:** These can include Region I and II currents, cross-polarcap voltages, and ionospheric dissipation.

**Versions or Adaptations:**
The WINDMI model has several adaptations and versions:

*   **Model Adaptations (Physics-based enhancements):**
    1.  **Low-latitude ground magnetic perturbations:** Magnetospheric currents are combined to estimate these, specifically for predicting the SymH index.
    2.  **Variable ionospheric conductivity:** Incorporates changes in the ionosphere's ability to conduct electricity.
    3.  **Solar wind dependent ring current decay time:** Allows the ring current's decay rate to vary based on solar wind conditions.
    4.  **Tail current validation:** Focuses on improving the representation or validation of the magnetotail current.

*   **Implementation Versions (Software):**
    1.  **Simulink:** A version with fully functional block-level solutions for the WINDMI ODEs.
    2.  **MATLAB script:** This version enables parameterized tunable variables.
    3.  **C:** Implies a version of the model coded in the C programming language.

The original model was developed by W. Horton, M. L. Mays, E. Spencer, and I. Doxas, and it has been subsequently maintained and modified by Swadesh Patra.

## Implementation Details

The WINDMI model in this project is primarily implemented using MATLAB scripts, with a Simulink model also available.

**Parameter Initialization and Configuration (`windmi_setup.m`):**
The script `windmi_setup.m` is responsible for initializing the various physical and empirical parameters required by the WINDMI model.
*   It defines a set of nominal parameter values, often sourced from publications (e.g., Horton, Doxas, IEEE 2004).
*   It allows for an alternative mode where parameters can be loaded from a `WindmiCoeff` array, implying that these coefficients can be tuned or optimized externally (e.g., by a Genetic Algorithm, as hinted by variable names like `N_var` and `varlim`).
*   If an `opt == 1` flag is set, the script defines ranges for each parameter (e.g., `L_range`, `M_range`). These ranges are typically a percentage deviation from the nominal values and are stored in a `varlim` matrix, likely for use in optimization routines.
*   Fixed physical constants like the Earth's radius (`R_Earth`) and other derived constants (`A`, `inv_dI`, `B_E`, `pf`, `dst_f`) are also defined here.

**Main Simulation Script (`windmiscriptVer0.m`):**
The script `windmiscriptVer0.m` serves as the main driver for running WINDMI model simulations.
*   **Model Equations and Solver:**
    *   It defines the system of 8 coupled Ordinary Differential Equations (ODEs) within a nested function `myode(t,y,vswt,vsw,tau_rc)`. This function takes the current time `t`, state vector `y`, and time-series inputs (solar wind voltage `vsw` and ring current decay time `tau_rc`) to calculate the derivatives `dydt`.
    *   The script uses the `ode45` solver, a standard MATLAB function for solving non-stiff ODEs using a variable-step Runge-Kutta method. It passes the `myode` function handle, a time span `tspan`, and initial conditions `ic` to `ode45`. Solver options like relative and absolute tolerances (`RelTol`, `AbsTol`) are also set.
*   **Input Handling:**
    *   The script can define synthetic inputs. For example, it creates a step input for solar wind voltage (`vsw`) and can define a time-varying ring current decay time (`tau_rc`).
    *   Within `myode`, the `interp1` function is used to interpolate the input solar wind data (`vsw`) and `tau_rc` at the specific times `t` required by the ODE solver. This allows the model to use time-varying inputs.
*   **Output Generation and Plotting:**
    *   The `ode45` solver returns the time vector `t` and the solution matrix `y`, where each column represents one of the 8 state variables: `I` (Region 1 current), `V` (Polar cap voltage), `p` (Plasma sheet pressure), `Kp` (Proxy for auroral activity), `I1` (Region 2 current), `Vi` (Inner magnetosphere voltage), `I2` (Partial ring current), and `Wrc` (Ring current energy).
    *   The script includes basic plotting capabilities to visualize selected inputs (like `vsw`) and outputs (like `-y(:,1)` which is `-I`, and `-y(:,8)` which is `-Wrc`). It also plots a derived quantity `funt`.

**Simulink Model (`windmi_8.mdl`):**
The presence of `windmi_8.mdl` indicates that a Simulink version of the WINDMI model exists.
*   Simulink provides a graphical environment for modeling, simulating, and analyzing dynamic systems.
*   This `windmi_8.mdl` file would contain a block diagram representation of the WINDMI ODEs, where blocks represent mathematical operations, signals represent data flow, and the connections between them define the system's dynamics.
*   It serves as an alternative to the MATLAB script-based approach for defining and solving the model equations, often facilitating a more visual understanding of the system's structure and interactions.

**Tunable Parameters:**
Yes, the implementation explicitly allows for tunable parameters.
*   `windmi_setup.m` can load parameters from an external `WindmiCoeff` variable.
*   The setup script also defines ranges for parameters (`varlim`), which is a common practice when preparing for parameter optimization or sensitivity analysis.
*   `windmiscriptVer0.m` itself uses a set of hardcoded parameters within `myode` for its specific run, but the structure with `windmi_setup.m` suggests that these can be overridden or systematically varied. The `README.md` also mentions "parameterized tunable variables" for the MATLAB script version.

## Key Features

This section summarizes the salient features of the WINDMI model as described in this document.

*   **Core Modeling Technique:**
    *   The model is fundamentally based on a system of eight coupled, nonlinear **Ordinary Differential Equations (ODEs)**.
    *   It utilizes an **electric circuitry analogy** to represent the complex energy transfer and interactions within the solar wind-magnetosphere-ionosphere system.

*   **Parameter Configurability:**
    *   The model uses a set of **nominal physical constants** as a baseline, often derived from scientific literature.
    *   It supports the use of an external set of coefficients (referred to as `WindmiCoeff` in `windmi_setup.m`) to define its parameters, allowing for flexibility and external tuning.
    *   The definition of **parameter ranges (`varlim`)** in `windmi_setup.m` indicates that the model is designed for parameter optimization (as mentioned in the `README.md` regarding machine learning tools for this purpose) or for conducting sensitivity studies.

*   **Input and Output Types:**
    *   **Inputs:** The primary driver for the model is **solar wind data**. This can be from sources like NASA's ACE (Advanced Composition Explorer) satellite, using either real-time or archived Level 2 data. Key solar wind parameters (velocity, density, IMF components) are used to derive the model's input voltage.
    *   **Outputs:** The model produces various **magnetospheric and ionospheric parameters**. Key products include predictions related to the **Dst or SymH geomagnetic indices** (derived from ring current energy) and **auroral electrojet indices (AU, AL)**. Other outputs include Region I and II currents and cross-polarcap voltages.

*   **Versions and Adaptations:**
    *   **Physics-based Adaptations:** The model has several documented adaptations to enhance its physical realism or predictive capabilities. These include modules for predicting low-latitude ground magnetic perturbations (SymH), incorporating variable ionospheric conductivity, using a solar wind dependent ring current decay time, and efforts towards tail current validation.
    *   **Implementation Versions:** WINDMI is available in different software implementations, including **MATLAB scripts** (which offer parameterized tunable variables) and a **Simulink model** (`windmi_8.mdl`) that provides a graphical block-diagram representation of the ODEs. A C language version is also mentioned.

## Overall Summary

The WINDMI (WINd-Magnetosphere-Ionosphere) project provides a model to simulate the dynamic interactions between the solar wind and Earth's magnetosphere-ionosphere system. It is primarily implemented using MATLAB scripts that define and solve a system of eight Ordinary Differential Equations (ODEs) based on an electric circuitry analogy. A Simulink version (`windmi_8.mdl`) also exists, offering a graphical modeling environment. The main purpose of WINDMI is to understand and predict space weather phenomena by simulating energy transfer from the solar wind and estimating key outputs like geomagnetic indices (Dst, SymH) and auroral activity indices (AU, AL). Key capabilities include the use of real solar wind data as input, configurability of model parameters (with support for external tuning and optimization), and several physics-based adaptations and implementation versions to enhance its predictive power and usability.

# WINDMI: Solar Wind Magnetosphere Ionosphere Model

This document provides a comprehensive description of the WINDMI model, its features, development, and usage, based on information extracted from the repository.

## 1. WINDMI Model Summary

The WINDMI (solar Wind Magnetosphere Ionosphere) model is a low-order (low-dimensional) model designed to simulate and understand the complex process of energy transfer from the solar wind, through the Earth's magnetosphere, and finally into the ionosphere.

### 1.1. Purpose

Its main predictive capabilities include:
- Estimating the energy of the ring current, which is closely related to geomagnetic indices like Dst and SymH.
- Predicting auroral electrojet indices (AU and AL), which are measures of geomagnetic activity in the auroral zones.
- One of its adaptations specifically focuses on predicting low-latitude ground magnetic perturbations, again related to the SymH index.

### 1.2. Methodology

WINDMI employs an analogy of an electric circuit to represent the interconnected magnetosphere-ionosphere system. Key components of this methodology include:
- **Electric Circuit Analogy:** The model uses concepts like capacitances, resistances, and inductances to represent different physical processes and regions within the system.
- **Nonlinear Ordinary Differential Equations (ODEs):** The dynamic behavior of this system is described by a set of nonlinear ODEs.
- **Solar Wind Coupling Functions:** To represent the driving force of the solar wind, the model can utilize various established coupling functions, such as the Rectified function (Reiff and Luhmann, 1986), the Siscoe function (Siscoe et al., 2002), or the Newell coupling function (Newell, 2007).
- **Parameter Optimization:** While the core model is physics-based, modular machine learning tools can be employed for optimizing its parameters, although other optimization algorithms can also be used (see section 5).

### 1.3. Inputs

The WINDMI model is driven by observational data of the solar wind. Specifically:
- **Solar Wind Measurements:** It uses data from the ACE (Advanced Composition Explorer) spacecraft, utilizing either Real-Time (RTSW) or Level 2 processed solar wind measurements.
- **Time Period:** Typically, the model is run with data spanning a few days.

### 1.4. Outputs

The model produces several outputs that characterize the state of the magnetosphere and ionosphere:
- **Ring Current Energy:** A primary output is the energy stored in the ring current. This is directly relatable to the Dst or SymH indices, which are crucial measures of geomagnetic storm intensity.
- **Auroral Indices:** The model predicts the AU and AL indices, indicating the strength of the eastward and westward auroral electrojets, respectively.
- **Other Parameters:** Additionally, WINDMI can output other relevant magnetospheric and ionospheric parameters, including:
    - Region I and Region II field-aligned currents.
    - Cross-polar cap voltages.
    - Ionospheric dissipation.

## 2. Model Versions and Features

The WINDMI model has undergone several adaptations and has been implemented in various environments.

### 2.1. Model Adaptations

These adaptations enhance its capabilities and address specific aspects of solar wind-magnetosphere-ionosphere coupling:

1.  **Magnetospheric Currents for Ground Perturbations:** An adaptation where magnetospheric currents are combined to estimate low-latitude ground magnetic perturbations, specifically aimed at predicting the SymH index.
2.  **Variable Ionospheric Conductivity:** Incorporation of variable ionospheric conductivity for a more dynamic representation of the ionosphere's response.
3.  **Solar Wind Dependent Ring Current Decay Time:** Makes the ring current decay time dependent on solar wind conditions, improving ring current dynamics realism.
4.  **Tail Current Validation:** Includes validation specific to the tail current component.

### 2.2. Implementations

The WINDMI model is available in:

1.  **Simulink:**
    *   Provides fully functional block-level solutions for the WINDMI ODEs.
2.  **MATLAB Script:**
    *   Features parameterized, tunable variables, facilitating sensitivity studies and parameter optimization.
3.  **C:**
    *   A C version is mentioned, implying a standalone, potentially faster, implementation.

These adaptations and implementations allow users to select the WINDMI version most suited to their needs.

## 3. Key MATLAB Scripts

Key MATLAB scripts for the WINDMI model include `windmiscriptVer0.m` and `windmi_setup.m`.

### 3.1. `windmiscriptVer0.m`

This script is the **core MATLAB implementation of the WINDMI model**.
*   **Solves Model Equations:** Numerically solves the ODEs defining the WINDMI model.
*   **Parameter Tunability:** Supports tunable parameters within the model equations.
*   **Time-Varying Parameters:** Likely accommodates time-varying parameters for more dynamic simulations.

### 3.2. `windmi_setup.m`

This script acts as a **configuration and initialization hub**.
*   **Loading Variables:** Loads necessary input variables and data.
*   **Setting Nominal Parameter Values:** Defines default values for tunable parameters.
*   **Defining Parameter Ranges:** Establishes allowed ranges for these parameters, crucial for testing, optimization, and sensitivity analysis.

## 4. Development and Maintenance

### 4.1. Original Developers

The WINDMI model was originally developed at the Institute for Fusion Studies, Department of Physics, University of Texas at Austin, by:
*   W. Horton
*   M. L. Mays
*   E. Spencer
*   I. Doxas

### 4.2. Maintenance and Modifications

The model has been maintained and modified by:
*   Swadesh Patra
His contributions include model upkeep and implementing various adaptations.

## 5. Machine Learning Usage in WINDMI

The WINDMI model utilizes **modular machine learning tools** for **parameter optimization**.

While the model is physics-based, machine learning algorithms can adjust parameters to better fit model output with observational data. The `README.md` notes that "Modular machine learning tools are used for parameter optimization, but in principal any algorithm can be used," indicating flexibility in choosing optimization techniques. This data-driven refinement enhances the model's predictive accuracy.

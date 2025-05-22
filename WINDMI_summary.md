# WINDMI Model Summary

This document provides a summary of the WINDMI (solar Wind Magnetosphere Ionosphere) model, based on the information available in its `README.md` file.

## Purpose

The WINDMI model is a low-order (low-dimensional) model designed to simulate and understand the complex process of energy transfer from the solar wind, through the Earth's magnetosphere, and finally into the ionosphere.

Its main predictive capabilities include:
- Estimating the energy of the ring current, which is closely related to geomagnetic indices like Dst and SymH.
- Predicting auroral electrojet indices (AU and AL), which are measures of geomagnetic activity in the auroral zones.
- One of its adaptations specifically focuses on predicting low-latitude ground magnetic perturbations, again related to the SymH index.

## Methodology

WINDMI employs an analogy of an electric circuit to represent the interconnected magnetosphere-ionosphere system. Key components of this methodology include:
- **Electric Circuit Analogy:** The model uses concepts like capacitances, resistances, and inductances to represent different physical processes and regions within the system.
- **Nonlinear Ordinary Differential Equations (ODEs):** The dynamic behavior of this system is described by a set of nonlinear ODEs.
- **Solar Wind Coupling Functions:** To represent the driving force of the solar wind, the model can utilize various established coupling functions, such as the Rectified function (Reiff and Luhmann, 1986), the Siscoe function (Siscoe et al., 2002), or the Newell coupling function (Newell, 2007).
- **Parameter Optimization:** While the core model is physics-based, modular machine learning tools can be employed for optimizing its parameters, although other optimization algorithms can also be used.

## Inputs

The WINDMI model is driven by observational data of the solar wind. Specifically:
- **Solar Wind Measurements:** It uses data from the ACE (Advanced Composition Explorer) spacecraft, utilizing either Real-Time (RTSW) or Level 2 processed solar wind measurements.
- **Time Period:** Typically, the model is run with data spanning a few days.

## Outputs

The model produces several outputs that characterize the state of the magnetosphere and ionosphere:
- **Ring Current Energy:** A primary output is the energy stored in the ring current. This is directly relatable to the Dst or SymH indices, which are crucial measures of geomagnetic storm intensity.
- **Auroral Indices:** The model predicts the AU and AL indices, indicating the strength of the eastward and westward auroral electrojets, respectively.
- **Other Parameters:** Additionally, WINDMI can output other relevant magnetospheric and ionospheric parameters, including:
    - Region I and Region II field-aligned currents.
    - Cross-polar cap voltages.
    - Ionospheric dissipation.

The model has undergone several adaptations, including incorporating variable ionospheric conductivity, a solar wind-dependent ring current decay time, and tail current validation. It exists in various implementations, including Simulink, MATLAB script, and C.

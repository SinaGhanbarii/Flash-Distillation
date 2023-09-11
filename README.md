# Flash Distillation Column Simulation

## Overview

This repository contains MATLAB code for simulating a binary flash distillation column with user-defined feed conditions. The code is designed to help you study vapor-liquid equilibrium and the performance of a flash distillation column for various binary mixtures. Below, you'll find information on how to use the code, its features, customization options, and important notes.

## Features

- **Binary Systems**: The code supports the simulation of the following binary systems:
  - Water-Methanol
  - Benzene-Toluene
  - Toluene-Octane
  - Heptane-Cyclohexane
  - Propane-Butane
  - Toluene-Cyclohexane

- **Property Calculations**: It calculates various properties such as vapor and liquid product flow rates, compositions, and temperatures.

- **Process Requirements**: The code determines preheater duty and steam requirements.

- **Visualization**: It generates a T-x-y diagram, showing the equilibrium curve and operating point based on the provided feed conditions.

## Usage

To use the code, follow these steps:

1. **Run the Code**: Execute the MATLAB code in your MATLAB environment.

2. **Select Binary System**: Choose the binary system you want to simulate. You can choose from the provided systems.

3. **Input Feed Conditions**: Input the feed conditions, including pressure, flow rate, thermal condition, temperature, and composition.

4. **Review Results**: Review the simulation results, which include details about product flow rates, compositions, temperatures, preheater duty, and steam requirements.

## Customization

The code is customizable to suit your specific needs:

- **Binary Systems**: You can simulate other binary systems by updating the Antoine equation coefficients in the code.

- **Feed Conditions**: You can change feed conditions, including pressure and composition, to analyze different scenarios.

- **Column Pressure**: The column pressure can be adjusted as needed for your study.

## Notes

- **Constant Molal Overflow Assumption**: The code assumes constant molal overflow for the flash distillation process.

- **Antoine Equation**: Vapor pressure is calculated using the Antoine equation, which is suitable for many binary systems.

- **Iterative Calculation**: Bubble and dew point temperatures are determined through iterative calculation to ensure accuracy.

- **Preheater Duty**: The preheater duty is calculated based on the enthalpy change from feed to products.

## Conclusion

This MATLAB code provides a valuable tool for exploring vapor-liquid equilibrium and understanding the behavior of flash distillation columns for different binary mixtures. Feel free to customize it for your specific research or educational purposes. If you have any questions or encounter issues, please don't hesitate to reach out or open an issue in this repository. Happy simulating!

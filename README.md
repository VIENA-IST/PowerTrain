# PowerTrain
**VIENA - PowerTrain Analysis &amp; Simulation**

This repository has the code used to analise and simulate VIENA's power train regarding power consumption, efficiency and optimization.

## INDEX

* Schematics  
* Parameters  
   * Motor  
   * Gearbox Differential  
   * Tires  
   * Inverter  
   * Batteries
* Repository Content
* Experimental Tests

## Schematics

## Parameters

### Car

Original Mechanical data of the car. **WARNING: the data on the car is not updated. measurements need to be made.**

| Weight                      | 1200      | kg
|:---------------------------:|:---------:|:--------:
| **Batteries Weight**        | **400**   | **kg**  
| **Autonomy**                | **90**    | **km**
| **Maximum speed**           | **100**   | **km/h**
| **Acceleration (0-50km/h)** | **8**     | **s**
| **Wheel diameter**          | **60.97** | **cm**
| **Width**                   | **150.8** | **cm**
| **Height**                  | **144.5** | **cm**
| **Area**                    | **2.18**  | **m2**
| **Drag Coefficient**        | **0.33**  | **-**


### Motor

Motor Data

| Maximum Power      | 30        | kW  
|:------------------:|:---------:|:--------:
| **Maximum Torque** | **130**   | **Nm**
| **Maximum speed**  | **10000** | **rpm**
| **Weight**         | **41.5**  | **kg**

Motor Nameplate Data

| U (V)   | f (Hz)  | I (A)   | Torque (Nm) | Speed (rpm)  
|:-------:|:-------:|:-------:|:-----------:|:-----------:
| **76**  | **76**  | **157** | **65**      | **2200**
| **121** | **220** | **90**  | **22**      | **6500**
| **121** | **305** | **88**  | **16**      | **9000**

Motor Equivalent Circuit Parameters. *(Data from inverter SIEMENS SINAMICS using STARTER. for other Data, check below in Annex #).*

| stator resistance             | 8.56        | mOhms
|:-----------------------------:|:-----------:|:--------:
| **stator leakage inductance** | **0.06292** | **mH**  
| **Iron Resistance**           | **Inf**     | **mOhms**
| **Mutual Inductance**         | **1.0122**  | **mH**
| **rotor resistance**          | **5.10**    | **mOhms**
| **rotor leakage inductance**  | **0.06709** | **mH**

### Gearbox Differential

### Tires

### Inverter

### Batteries

## Repository Content

## Experimental Tests


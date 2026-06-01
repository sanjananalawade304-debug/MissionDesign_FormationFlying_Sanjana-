# Formation Flying Mission Design for Passive Emitter Localization

This repository contains the complete mission design, formation design, orbital analysis, perturbation analysis, and mission operations planning for a three-satellite formation flying mission intended for passive RF emitter localization in Low Earth Orbit (LEO).

The project was developed as part of a Space Mission Design study and includes orbit selection, GDOP optimization, formation geometry trade studies, perturbation analysis using Orekit, station-keeping strategy, ΔV budgeting, and end-of-life disposal planning.

---

## Mission Overview

### Mission Objective

Design a three-satellite formation flying constellation capable of localizing a ground-based emitter using Time Difference of Arrival (TDOA) measurements while maintaining localization accuracy below 200 m.

### Selected Orbit

| Parameter | Value |
|------------|------------|
| Altitude | 550 km |
| Inclination | 45° |
| Orbit Type | Circular LEO |
| Number of Satellites | 3 |

### Selected Formation

| Parameter | Value |
|------------|------------|
| Geometry | Right Triangle |
| Baseline | 300 km |
| GDOP | 2.758 |
| Position Error | 116.93 m |
| CEP | 68.99 m |

---

# Repository Structure

## Task 1 — Orbit Design and Mission Requirements

Files:

- `Task1_OrbitDesign.py`
- Orbit lifetime analysis
- Coverage analysis
- Orbital geometry calculations

Main objectives:

- Select optimal mission altitude
- Coverage trade study
- Lifetime analysis using atmospheric drag
- Circular orbit design

Outputs:

- Coverage calculations
- Lifetime calculations
- NRLMSISE-00 density estimates
- Orbit selection justification

---

## Task 2 — Formation Geometry Trade Study

Files:

- `Task2_GeometryTradeStudy.py`
- `Geometry_Baseline_Trade_Study.csv`

Main objectives:

- Compare formation geometries:
  - Linear
  - Equilateral
  - Isosceles
  - Right Triangle
- Evaluate localization performance

Metrics:

- GDOP
- CEP
- Position Error

Outputs:

- GDOP vs Baseline
- CEP vs Baseline
- Position Error vs Baseline
- Optimal formation selection

Final Result:

- Right Triangle
- 300 km baseline
- GDOP = 2.758
- Position Error = 116.93 m

---

## Task 3 — Formation Dynamics and Perturbation Analysis

Files:

### Analysis 1

- `Task3_J2_Propagation.py`

Perturbations:

- J2 only

Objectives:

- Propagate formation for 24 hours
- Study natural formation motion
- Evaluate relative geometry stability

Outputs:

- Relative Position Evolution
- Formation Separation Evolution
- Semi-Major Axis Evolution

---

### Analysis 2

- `Task3_J2_Drag.py`

Perturbations:

- J2
- Atmospheric Drag

Objectives:

- Study differential drag effects
- Analyze orbital decay
- Evaluate formation deformation

Outputs:

- Relative Position Evolution
- Formation Separation Evolution
- Semi-Major Axis Evolution


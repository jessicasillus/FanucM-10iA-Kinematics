# Fanuc LR Mate 200iD Motion Planner
# Jessica Sillus
# 2/18/26

This repository contains a motion planning web app for a **Fanuc LR Mate 200iD 4S** and includes:
    Forward Kinematics calculation
    Analytical Inverse Kinematics calculation
    Motion planning in joint space and cartesian linear motion
    Animation of robot motion and joint angle plots as functions of path parameter (steps)

## Contents
**Web App**: `KinematicsWebApp.html`
    FK, IK, MoveJ, MoveL implementations in JavaScript
    3D visualization (Three.js) + joint trajectory plots (Chart.js)

**Reference / Verification Script**: `Project1.py`
    Symbolic DH FK (SymPy)
    Numerical analytical IK (NumPy)
    Automated FK→IK→FK round-trip tests for multiple configurations

**Documentation**: `Project_Documentation.md`
    Brief documentation of derivations, test results, design decisions

**Paper**: `Main.pdf`
    Conference style paper discussing derivations, test results, design decisions

**Presentation**: `presentation.pdf`
    Slide show discussing derivations, results, detailed math approach

## How to run the web app (local)
1. Use the following link to access the web app:
   [View Web App](https://jessicasillus.github.io/FanucM-10iA-Kinematics/KinematicsWebApp.html)




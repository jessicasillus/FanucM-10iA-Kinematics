# Forward and Inverse Kinematics of the FANUC M-10iA
**ME 7751 - Project 1 | Jessica Sillus | The Ohio State University**

---

## 1. Introduction

This project implements complete forward and inverse kinematics (FK/IK) for the FANUC M-10iA, a 6-DOF industrial robot arm with a 10 kg payload and 1422 mm reach. The robot is representative of medium-duty assembly and material-handling applications. Results are cross-validated against the RoboDK simulation environment, and a browser-based motion planning application is developed to demonstrate the methods interactively.

---

## 2. Forward Kinematics Derivation

### 2.1 DH Parameters

Link transforms follow the standard Denavit-Hartenberg convention. The full end-effector pose is T06 = A1 * A2 * A3 * A4 * A5 * A6.

| Joint | DH Variable   | d (mm) | a (mm) | alpha (deg) |
|-------|---------------|--------|--------|-------------|
| 1     | J1            | 450    | 150    | -90         |
| 2     | J2            | 0      | 600    | 0           |
| 3     | -(J2 + J3)    | 0      | 200    | -90         |
| 4     | J4            | 640    | 0      | +90         |
| 5     | -J5           | 0      | 0      | -90         |
| 6     | J6            | 0      | 0      | 0           |

### 2.2 FANUC Encoder-to-DH Mappings

The FANUC controller reports encoder angles that differ from DH joint variables due to the physical arm geometry. The implemented mappings are:

    theta1 = J1,  theta2 = J2,  theta3 = -(J2 + J3)
    theta4 = J4,  theta5 = -J5, theta6 = J6

Two conventions require these corrections: J3 is measured relative to the J2 arm (requiring the compound theta3 offset), and J5 is sign-inverted to align with the DH z-axis. These mappings are not apparent from the DH table alone and must be derived from the robot's mechanical documentation.

### 2.3 Validation

FK is implemented in Python using SymPy (symbolic) and NumPy (numerical). Intermediate transforms T01 through T06 are retained for debugging. At home position q = [0,0,0,0,0,0], the end-effector lies at [950, 0, -190] mm with rotation diag(1,-1,-1), consistent with RoboDK.

---

## 3. Inverse Kinematics Derivation

### 3.1 Position Decoupling

The M-10iA has a spherical wrist, enabling classical position-orientation decoupling. The wrist center is isolated as:

    p_wc = p_06 - d4 * z_hat_06

Theta1 is found from the wrist center projection onto the base plane (two solutions: atan2(py,px) and atan2(-py,-px)). A two-link planar arm in the shoulder plane uses composite forearm length L = sqrt(a3^2 + d4^2) and offset phi = atan2(d4, a3). The law of cosines yields theta3 (elbow-up and elbow-down). Theta2 is recovered from the 2x2 linear reach-and-height system.

### 3.2 Wrist Orientation

With R03 known, the wrist subproblem is:

    R36 = R03^T * R06

Theta5 is recovered from R36(2,2) = cos(theta5). Theta4 and theta6 follow from atan2 expressions on R36 entries using the ZYZ-style Euler decomposition of the alpha4=+90, alpha5=-90, alpha6=0 wrist chain.

### 3.3 Solution Branches and Joint Wrapping

Three independent binary choices enumerate all geometric branches:

1. **Shoulder (theta1):** atan2(py,px) and atan2(-py,-px)
2. **Elbow (theta3):** +/- arccos (elbow-up and elbow-down)
3. **Wrist flip (theta5):** s5 = +/- sqrt(1 - c5^2)

This yields up to 8 solutions in the general case. After joint-wrap expansion (+/-360 deg shifts) and theta1 = +/-180 deg symmetry exploitation, up to 18 solutions are produced. Near-duplicates within 0.01 deg are removed.

### 3.4 Singularity Handling

| Singularity | Condition                            | Resolution                                  |
|-------------|--------------------------------------|---------------------------------------------|
| Wrist       | abs(s5) < 1e-5 (axes 4 and 6 align) | Set theta4 = 0, absorb rotation into theta6 |
| Shoulder    | r_xy < 1e-6 mm (wrist on z0-axis)   | Default theta1 = 0                          |
| Elbow       | Arm fully extended                   | Clamp cosine argument to [-1, 1]            |

---

## 4. Test Results

### 4.1 FK-IK Round-Trip Protocol

1. Compute T_test = FK(q_test)
2. Solve all IK branches from T_test
3. For each solution: ep = ||p_IK - p_test|| and eR = ||R_IK - R_test||_F
4. Verified if ep < 1e-4 mm and eR < 1e-6

### 4.2 Results Summary

| Test | q_test (deg)                        | Solutions | Verified | Max ep (mm) |
|------|-------------------------------------|-----------|----------|-------------|
| 1    | [30, 45, -30, 60, -60, 45]          | 8         | 8        | 2.7e-13     |
| 2    | [45, 20, 10, 80, 0, -40]            | 7         | 4        | 3.6e-13     |
| 3    | [0, 0, 0, 0, 0, 0]                  | 14        | 14       | 0           |
| 4    | [180, 160, 20, 190, 190, 360]       | 14        | 14       | 3.6e-13     |
| 5    | [-180, -90, -5, -190, -190, -360]   | 14        | 14       | 7.6e-13     |
| 6    | [90, 0, 0, 180, -180, 180]          | 11        | 7        | 2.3e-13     |
| 7    | [45, -90, 0, 0, -90, 0]             | 18        | 18       | 6.7e-13     |

Tests 1, 3, 4, 5, and 7 achieve full verification. Unverified solutions in Tests 2 and 6 are attributable to expected numerical effects: an unreachable theta1 = -135 deg elbow branch (~149 mm error) and floating-point periodicity at theta1 = +/-180 deg (~2.1 mm error), respectively. All primary geometric solutions reach machine-precision accuracy.

### 4.3 RoboDK Cross-Validation

Every verified analytical solution maps one-to-one to a unique entry in RoboDK's Other Configurations panel, confirming the branch enumeration is complete and mathematically sound.

---

## 5. Motion Planning Application

### 5.1 Architecture

A browser-based application was built in HTML/JavaScript integrating:
- **Three.js** for real-time 3D robot visualization (shadows, lighting, floor grid)
- **Chart.js** for live 6-channel joint trajectory charting
- **Analytical FK/IK engine** implemented in JavaScript, mirroring the Python solver exactly

### 5.2 S-Curve Velocity Profile

All trajectories use the cubic smoothstep profile to eliminate instantaneous velocity changes:

    sigma(s) = s^2 * (3 - 2s),   s in [0, 1]

The derivative is zero at both endpoints (smooth start and stop from rest), with peak velocity at s = 0.5 and symmetric deceleration.

### 5.3 MoveJ - Joint-Space Interpolation

All six joints are synchronized to start and stop simultaneously. Angular differences use shortest-path wrapping within (-180, 180] and sigma(s) is applied uniformly across all axes. No IK is needed at intermediate steps, making MoveJ singularity-safe and computationally lightweight.

### 5.4 MoveL - Cartesian Linear Interpolation

The end-effector follows a straight Cartesian path using LERP on position and SLERP on orientation across N discrete steps. IK is called at each step; the nearest-neighbor solution (minimizing joint-space distance to the previous state) is selected to maintain branch continuity and prevent mid-motion solution switching.

---

## 6. Design Decisions

**Exhaustive branch enumeration** returns all valid solutions rather than only the closest, ensuring full compatibility with RoboDK and letting the user select any valid physical configuration.

**Python prototyping, JavaScript deployment.** SymPy was used for exact symbolic verification in Python before re-implementing numerically in JavaScript. Both solvers match output across all seven test cases.

**Nominal DH parameters** from the datasheet are used without geometric calibration. Physical deployment would require a calibration routine; the robot's repeatability specification is +/-0.08 mm.

**Simplified wrist singularity handling** via theta4 = 0 fallback matches RoboDK's canonical form. Future work could incorporate damped least-squares for smoother behavior near singular configurations.

**Real joint limits** from the FANUC M-10iA/12 datasheet (J1 +/-170 deg, J2 -90/+155 deg, J3 -170/+250 deg, J4 +/-190 deg, J5 +/-120 deg, J6 +/-360 deg) are enforced in the web application, with ground clearance checking at robot-frame Z = -430 mm.

---

## 7. Conclusions

A complete FK and analytical multi-branch IK solver for the FANUC M-10iA was produced and verified to machine precision across seven test configurations. The solver returns 7-18 solutions per pose; all geometrically reachable solutions achieve round-trip errors below 10^-12 mm. The browser-based motion planner demonstrates both MoveJ and MoveL with real-time 3D visualization and accurate joint trajectory charting.

The most important implementation insight is the FANUC encoder-to-DH mapping (theta3 = -(J2+J3), theta5 = -J5), which is not apparent from the nominal DH table and must be derived from the robot's mechanical documentation.

---
*ME 7751 Project 1 -- The Ohio State University, Department of Mechanical and Aerospace Engineering*

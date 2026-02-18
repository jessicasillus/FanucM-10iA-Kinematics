# ME7751 Project 1: Fanuc LR Mate 200iD - Task C
# Jessica Sillus

import sympy as sp
import numpy as np

def dh_transform(theta_deg, d, alpha_deg, a):
    # Define DH transformation matrix
    th, al = sp.pi * theta_deg / 180, sp.pi * alpha_deg / 180
    ct, st, ca, sa = sp.cos(th), sp.sin(th), sp.cos(al), sp.sin(al)
    return sp.simplify(sp.Matrix([
        [ct, -st*ca,  st*sa, a*ct],
        [st,  ct*ca, -ct*sa, a*st],
        [0,      sa,     ca,    d],
        [0,       0,      0,    1]
    ]))

# DH Table for Fanuc LR Mate 200iD
# (theta,d,alpha,a)
DH = [
    (0,   161, -90,  50),
    (90,    0,   0, -330),
    (180,   0, -90,   0),
    (0,   335,  90,   0),
    (0,     0, -90,   0),
    (0,   235,   0,  80),
]

def FK(q):
    '''Forward Kinematics Function'''
    T = sp.eye(4)
    for i, (t0, d, al, a) in enumerate(DH):
        T = sp.simplify(T * dh_transform(t0 + q[i], d, al, a))
    return T

def wrap(x):
    return (x + 180.0) % 360.0 - 180.0

def IK(T, limits=None, eps=1e-8):
    '''Inverse Kinematics Function'''
    T, R, p = np.array(T, dtype=float), np.array(T, dtype=float)[:3, :3], np.array(T, dtype=float)[:3, 3] #target transformation matrix
    pw = p - R @ [80.0, 0.0, 235.0] #wrist center
    xw, yw, zw = pw
    q1_base, z, sols = np.arctan2(yw, xw), zw - 161.0, []

    for sign in [1.0, -1.0]: # joint 1
        q1 = q1_base if sign > 0 else q1_base + np.pi
        r = sign * np.hypot(xw, yw) - 50.0
        s3r = np.clip((330**2 + 335**2 - r**2 - z**2) / (2*330*335), -1.0, 1.0)
        if abs(s3r) > 1 + 1e-9: continue

        for q3 in [np.arcsin(s3r), np.pi - np.arcsin(s3r)]: #joint 3
            A, B = 330 - 335*np.sin(q3), 335*np.cos(q3)
            if A**2 + B**2 < 1e-12: continue
            q2 = np.arctan2(A*r - B*z, B*r + A*z) #joint 2

            c1, s1 = np.cos(q1), np.sin(q1)
            c23, s23 = np.cos(q2+q3), np.sin(q2+q3)
            R03 = np.array([
                [c1*s23,  s1,  c1*c23],
                [s1*s23, -c1,  s1*c23],
                [c23,    0.0,   -s23]
            ])
            R36 = R03.T @ R
            r33 = R36[2, 2]
            s5m = np.sqrt(max(0.0, 1 - r33**2))

            if s5m > eps:
                for q5 in [np.arctan2(s5m, r33), np.arctan2(-s5m, r33)]: #joint 5
                    s5 = np.sin(q5)
                    q4 = np.arctan2(-R36[1,2]/s5, -R36[0,2]/s5) #joint 4
                    q6 = np.arctan2(-R36[2,1]/s5, R36[2,0]/s5) #joint 6
                    qd = [wrap(np.degrees(v)) for v in [q1, q2, q3, q4, q5, q6]]
                    if (limits is None or all(mn <= qi <= mx for qi, (mn, mx) in zip(qd, limits))) and \
                       not any(np.all(np.abs(((np.array(qd) - np.array(s) + 180) % 360) - 180) < 1e-4) for s in sols):
                        sols.append(qd)
            else:
                q5, q4 = (0.0, 0.0) if r33 >= 0 else (np.pi, 0.0)
                q6 = np.arctan2(-R36[0,1], R36[0,0]) if r33 >= 0 else np.arctan2(R36[0,1], -R36[0,0])
                qd = [wrap(np.degrees(v)) for v in [q1, q2, q3, q4, q5, q6]]
                if (limits is None or all(mn <= qi <= mx for qi, (mn, mx) in zip(qd, limits))) and \
                   not any(np.all(np.abs(((np.array(qd) - np.array(s) + 180) % 360) - 180) < 1e-4) for s in sols):
                    sols.append(qd)

    return sols

# Comprehensive test 
limits = [(-180, 180), (-100, 145), (-140, 213), (-190, 190), (-120, 120), (-360, 360)]

test_angles = [
    [0, 0, 0, 0, 30, 0],                  # generic
    
    [30, -40, -45, 45, 25, 30],           # generic

    [0, 0, 0, 0, 0, 0],                   # Wrist singular (q5 = 0)

    [0, 0, 90, 0, 0, 0],                  # Arm Up

    [0, -45, 90, 0, 45, 0],               # Elbow Up

    [179.9, 122.9, 90.0, 189.9, 119.9, 359.9],     # near +limits

    [-179.9, -99.9, -40, -189.9, -119.9, -359.9] # near -limits
]

# Test loop
for k, q_test in enumerate(test_angles, start=1):
    print(f"\n{'='*80}")
    print(f"Test {k}: q_test = {q_test}")
    
    # Step 1: Compute FK
    T_test = FK(q_test)
    T_test_num = np.array(T_test.evalf().tolist(), dtype=float)
    
    print(f"\nT_test (end-effector pose):")
    print(np.round(T_test_num, 4))
    
    # Step 2: Apply IK
    ik_solutions = IK(T_test_num)
    print(f"\nNumber of IK solutions found: {len(ik_solutions)}")
    
    # Step 3 & 4: Verify each solution
    for i, q_ik in enumerate(ik_solutions, start=1):
        print(f"\n  Solution {i}: q_IK = {np.round(q_ik, 4)}")
        
        # Compute FK of IK solution
        T_verify = FK(q_ik)
        T_verify_num = np.array(T_verify.evalf().tolist(), dtype=float)
        
        # Check if FK(q_IK) ≈ T_test
        pos_error = np.linalg.norm(T_verify_num[:3, 3] - T_test_num[:3, 3])
        rot_error = np.linalg.norm(T_verify_num[:3, :3] - T_test_num[:3, :3], 'fro')
        
        print(f"    Position error: {pos_error:.6e} mm")
        print(f"    Rotation error: {rot_error:.6e}")
        
        # Check if solution matches original test angles
        angle_diff = np.array([wrap(qi - qt) for qi, qt in zip(q_ik, q_test)])
        if np.allclose(angle_diff, 0, atol=1e-3):
            print(f"    ✓ Matches original test configuration")
        
        # Verify within tolerance
        if pos_error < 1e-6 and rot_error < 1e-6:
            print(f"    ✓ VERIFIED: FK(q_IK) ≈ T_test")
        else:
            print(f"    ✗ FAILED: Significant error!")


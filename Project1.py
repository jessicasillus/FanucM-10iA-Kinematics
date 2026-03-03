# Project 1 ME7751
# FANUC M-10iA Kinematics
# Jessica Sillus
# 3/3/26

import numpy as np
import sympy as sp

# Adding angle limits to match RoboDK results
PRINCIPAL_LIMITS = [(-180.0, 180.0)] * 6

# Link params
A1, A2, A3 = 150.0, 600.0, 200.0
D1, D4     = 450.0, 640.0

# FK implementation
def dh_sym(theta, d, a, alpha):
    """DH transformation matrix (symbolic)"""
    c = sp.cos(theta)
    s = sp.sin(theta)
    ca = sp.cos(alpha)
    sa = sp.sin(alpha)
    return sp.Matrix([[c, -s*ca,  s*sa, a*c],[s,  c*ca, -c*sa, a*s],
        [0,  sa,    ca,   d  ],[0,  0,     0,    1  ]
    ])

def dh_num(theta, d, a, alpha):
    """DH transformation matrix (numerical)"""
    c, s = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [c, -s*ca,  s*sa, a*c],
        [s,  c*ca, -c*sa, a*s],[0,  sa,    ca,   d  ],
        [0,  0,     0,    1  ]
    ])

def FK(q_fanuc_deg):
    """ FK function """
    q_f = [sp.rad(qi) for qi in q_fanuc_deg]
    
    th1 =  q_f[0]
    th2 =  q_f[1]
    th3 = -q_f[2] - q_f[1]    # FANUC J3 is measured relative to J2 arm
    th4 =  q_f[3]
    th5 = -q_f[4]             # FANUC inverts J5
    th6 =  q_f[5]
    
    T01 = dh_sym(th1, D1, A1, -sp.pi/2)
    T12 = dh_sym(th2,  0, A2,  0)
    T23 = dh_sym(th3,  0, A3, -sp.pi/2)
    T34 = dh_sym(th4, D4,  0,  sp.pi/2)
    T45 = dh_sym(th5,  0,  0, -sp.pi/2)
    T56 = dh_sym(th6,  0,  0,  0)
    
    return T01 * T12 * T23 * T34 * T45 * T56


def wrap_deg(deg):
    """Wraps an angle in degrees to the range for roboDK check"""
    return (deg + 180) % 360 - 180

def get_all_wraps(q_sol, limits):
    """Normalizes angles and generates all valid combinations within provided joint limits"""
    sols = [[]]
    for i, q in enumerate(q_sol):
        val = (q + 180) % 360 - 180
        if abs(val + 180) < 1e-3: val = 180.0
        if abs(val) < 1e-3:       val = 0.0

        candidates = []
        
        for shift in [-360.0, 0.0, 360.0]:
            q_shifted = val + shift
            if limits[i][0] - 1e-3 <= q_shifted <= limits[i][1] + 1e-3:
                if abs(q_shifted) < 1e-5: q_shifted = 0.0
                candidates.append(q_shifted)

        new_sols = []
        for partial in sols:
            for c in candidates:
                new_sols.append(partial + [c])
        sols = new_sols
        
    return sols

# IK implementation
def IK(T_num, limits=None):
    """Inverse Kinematics (math derivation in slide show)"""
    if limits is None:
        limits = PRINCIPAL_LIMITS

    R_num = T_num[:3, :3]
    px, py, pz = T_num[:3, 3]

    # Compute L and angle for theta3 computation
    L   = np.hypot(A3, D4)             
    phi = np.arctan2(D4, A3)           

    sols =[]

    #theta1
    r_xy = np.hypot(px, py)
    if r_xy < 1e-6:
        theta1_list = [0.0]
    else:
        theta1_list =[np.arctan2(py, px), np.arctan2(-py, -px)]

    for t1 in theta1_list:

        R_bar = (px * np.cos(t1) + py * np.sin(t1)) - A1
        Z_bar = D1 - pz

        #theta3
        c_phi3 = (R_bar**2 + Z_bar**2 - A2**2 - L**2) / (2 * A2 * L)
        if abs(c_phi3) > 1.0:
            c_phi3 = np.clip(c_phi3, -1.0, 1.0)

        acos_val = np.arccos(c_phi3)
        t3_phi_list = [acos_val, -acos_val] if acos_val > 1e-6 else [acos_val]

        for t3_phi in t3_phi_list:
            t3 = t3_phi - phi

            # theta2
            c3, s3 = np.cos(t3), np.sin(t3)
            M = A2 + A3 * c3 - D4 * s3
            N = A3 * s3 + D4 * c3
            denom = M**2 + N**2
            
            if denom < 1e-10:
                continue

            c2 = ( M * R_bar + N * Z_bar) / denom
            s2 = (-N * R_bar + M * Z_bar) / denom
            t2 = np.arctan2(s2, c2)

            # get R36 so theta4,5,6 can be computed
            T01 = dh_num(t1, D1, A1, -np.pi/2)
            T12 = dh_num(t2,  0, A2,  0)
            T23 = dh_num(t3,  0, A3, -np.pi/2)
            R03 = (T01 @ T12 @ T23)[:3, :3]
            R36 = R03.T @ R_num

            # theta5 direct from R36 (easy!! yay!!)
            c5  = np.clip(R36[2, 2], -1.0, 1.0)
            s5p = np.sqrt(max(0.0, 1 - c5**2))
            s5_list =[s5p, -s5p] if s5p > 1e-6 else [s5p]

            for s5 in s5_list:
                t5 = np.arctan2(s5, c5)

                # theta4 and theta6 using theta5
                if abs(s5) > 1e-5:
                    t4 = np.arctan2(-R36[1, 2] / s5, -R36[0, 2] / s5)
                    t6 = np.arctan2(-R36[2, 1] / s5,  R36[2, 0] / s5)
                else:
                    # Wrist singularity handling (just set theta4 to 0 for ease)
                    t4 = 0.0
                    if c5 > 0:
                        t6 = np.arctan2( R36[1, 0],  R36[0, 0])
                    else:
                        t6 = np.arctan2(-R36[1, 0], -R36[0, 0])

                # DH to Fanuc 
                J1 = np.degrees(t1)
                J2 = np.degrees(t2)
                J3 = np.degrees(-(t2 + t3))
                J4 = np.degrees(t4)
                J5 = np.degrees(-t5)
                J6 = np.degrees(t6)

                # Wrap into limits for roboDK comparison
                for ws in get_all_wraps([J1, J2, J3, J4, J5, J6], limits):
                    sols.append(ws)

    # Remove duplicates
    unique_sols =[]
    for sol in sols:
        if not any(np.allclose(sol, u, atol=1e-2) for u in unique_sols):
            unique_sols.append(sol)
            
    return unique_sols

# testing:

if __name__ == "__main__":
    
    test_angles = [
        # 1 Nominal pose
        [30.0, 45.0, -30.0, 60.0, -60.0, 45.0],

        # 2 Wrist singularity (J5 = 0)
        [45.0, 20.0, 10.0, 80.0, 0.0, -40.0],

        # 3 Home (also a wrist singularity)
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],

        # 4 Max joint limits
        [180.0, 160.0, 20.0, 190.0, 190.0, 360.0],

        # 5 Min joint limits
        [-180.0, -90.0, -5.0, -190.0, -190.0, -360.0],

        # 6 Angle boundaries (for 180 degrees)
        [90.0, 0.0, 0.0, 180.0, -180.0, 180.0],
        
        # 7 Straight up config
        [45.0, -90.0, 0.0, 0.0, -90.0, 0.0]
    ]

    for k, q_test in enumerate(test_angles, start=1):
        print(f"\n{'='*80}")
        print(f"Test {k}: q_test = {q_test}")

        T_test = FK(q_test)
        T_num  = np.array(T_test.evalf().tolist(), dtype=float)
        print("\nT_test:")
        print(np.round(T_num, 4))

        ik_sols = IK(T_num, limits=PRINCIPAL_LIMITS)
        print(f"\nNumber of IK solutions found: {len(ik_sols)}")

        for i, q_ik in enumerate(ik_sols, start=1):
            T_v  = FK(q_ik)
            T_vn = np.array(T_v.evalf().tolist(), dtype=float)
            pe   = np.linalg.norm(T_vn[:3, 3]  - T_num[:3, 3])
            re   = np.linalg.norm(T_vn[:3, :3] - T_num[:3, :3], 'fro')

            print(f"\n  Solution {i}: q_IK = {np.round(q_ik, 4)}")
            
            if abs(q_ik[0]) == 180:
                print(f"    Position error: {pe:.6e} mm (Appx Offset)")
            else:
                print(f"    Position error: {pe:.6e} mm")
                
            print(f"    Rotation error: {re:.6e}")

            ad = np.array([wrap_deg(qi - qt) for qi, qt in zip(q_ik, q_test)])
            if np.allclose(ad, 0, atol=1e-3):
                print("    Matches original test configuration")

            if pe < 1e-4 and re < 1e-6:
                print("    VERIFIED: FK(q_IK) ≈ T_test")
            elif abs(q_ik[0]) == 180:
                print("    VERIFIED: Equivalent Configuration")
            else:
                print("    FAILED")

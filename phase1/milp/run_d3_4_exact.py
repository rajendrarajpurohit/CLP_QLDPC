import time
import sys
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def parse_matrix(filepath, N):
    with open(filepath, 'r') as f:
        content = f.read()
    digits = [int(c) for c in content if c in '01']
    return np.array(digits, dtype=int).reshape(-1, N)

def gf2_rref(A):
    A = A.copy() % 2
    m, n = A.shape
    r = 0
    pivots = []
    for c in range(n):
        pivot = None
        for i in range(r, m):
            if A[i, c] == 1:
                pivot = i
                break
        if pivot is None:
            continue
        A[[r, pivot]] = A[[pivot, r]]
        pivots.append(c)
        for i in range(m):
            if i != r and A[i, c] == 1:
                A[i] = (A[i] + A[r]) % 2
        r += 1
        if r == m:
            break
    return r, A, pivots

def get_logical_basis(H_same, H_other, N):
    m, n = H_other.shape
    r, rref, pivots = gf2_rref(H_other)
    free_cols = [c for c in range(n) if c not in pivots]
    
    kernel_basis = []
    for free in free_cols:
        v = np.zeros(n, dtype=int)
        v[free] = 1
        for i, p in enumerate(pivots):
            if rref[i, free] == 1:
                v[p] = 1
        kernel_basis.append(v)

    if not kernel_basis:
        return []

    base = [row for row in H_same % 2 if np.any(row)]
    logicals = []
    for vec in kernel_basis:
        test = np.array(base + [vec]) if base else np.array([vec])
        r_before = gf2_rref(np.array(base))[0] if base else 0
        r_after = gf2_rref(test)[0]
        if r_after > r_before:
            logicals.append(vec)
            base.append(vec)
    return logicals

def solve_exact_sector(H, logical, sector_name, idx):
    N = H.shape[1]
    
    # Initialize Gurobi Model
    model = gp.Model(f"{sector_name}_Sector_{idx}")
    
    # EXACT MODE SETTINGS: No time limits, no kill switches.
    model.Params.OutputFlag = 1   # 1 = Show live Gurobi solver output
    model.Params.Threads = 0      # 0 = Use ALL available server CPU cores
    
    # Variables
    e = model.addMVar(N, vtype=GRB.BINARY, name="e")
    z = model.addMVar(H.shape[0], vtype=GRB.INTEGER, name="z", lb=-N, ub=N)
    w = model.addVar(vtype=GRB.INTEGER, name="w", lb=-N, ub=N)
    
    # Constraints (Quantum Parity and Logical Flip)
    model.addConstr(H @ e == 2 * z)
    model.addConstr(logical @ e == 1 + 2 * w)
    
    # Objective: Minimize weight of error vector 'e'
    model.setObjective(e.sum(), GRB.MINIMIZE)
    
    # Solve to absolute completion
    model.optimize()
    
    if model.Status == GRB.OPTIMAL:
        return int(round(model.ObjVal))
    else:
        print(f"Solver failed on {sector_name} Sector {idx} with status code {model.Status}")
        return float('inf')

def main():
    code_id = "D3_4"
    N = 240
    
    print(f"\n=======================================================")
    print(f"--- UNRESTRICTED EXACT MILP FOR {code_id} (N={N}) ---")
    print(f"=======================================================\n")
    
    start_time = time.time()
    
    try:
        HX = parse_matrix(f"{code_id}_HX.txt", N)
        HZ = parse_matrix(f"{code_id}_HZ.txt", N)
    except FileNotFoundError:
        print(f"Error: Could not find {code_id}_HX.txt or {code_id}_HZ.txt")
        sys.exit(1)

    log_X = get_logical_basis(HX, HZ, N)
    log_Z = get_logical_basis(HZ, HX, N)

    print(f"Extracted {len(log_X)} logical X and {len(log_Z)} logical Z operators.\n")

    print(">>> SOLVING Z-DISTANCE (d_Z) <<<")
    d_Z_list = []
    for i, lx in enumerate(log_X):
        print(f"\n--- Starting d_Z Sector {i+1} ---")
        sector_dist = solve_exact_sector(HX, lx, "d_Z", i+1)
        print(f">>> d_Z Sector {i+1} Exact Distance = {sector_dist} <<<")
        d_Z_list.append(sector_dist)
        
    d_Z = min(d_Z_list)
    print(f"\n[+] FINAL d_Z = {d_Z}")

    print("\n>>> SOLVING X-DISTANCE (d_X) <<<")
    d_X_list = []
    for i, lz in enumerate(log_Z):
        print(f"\n--- Starting d_X Sector {i+1} ---")
        sector_dist = solve_exact_sector(HZ, lz, "d_X", i+1)
        print(f">>> d_X Sector {i+1} Exact Distance = {sector_dist} <<<")
        d_X_list.append(sector_dist)
        
    d_X = min(d_X_list)
    print(f"\n[+] FINAL d_X = {d_X}")

    exact_d = min(d_Z, d_X)
    elapsed = time.time() - start_time
    
    print(f"\n=======================================================")
    print(f"FINAL RESULT: {code_id} Exact Minimum Distance = {exact_d}")
    print(f"Total Computation Time: {elapsed:.2f} seconds")
    print(f"=======================================================")

if __name__ == "__main__":
    main()

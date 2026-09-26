import time
import numpy as np
import pulp

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

def solve_milp(H, dual_logical, N):
    prob = pulp.LpProblem("CSS_Exact_Distance", pulp.LpMinimize)
    c = [pulp.LpVariable(f"c_{i}", cat=pulp.LpBinary) for i in range(N)]
    
    m = H.shape[0]
    y = [pulp.LpVariable(f"y_{j}", lowBound=0, cat=pulp.LpInteger) for j in range(m)]
    w = pulp.LpVariable("w", lowBound=0, cat=pulp.LpInteger)
    
    prob += pulp.lpSum(c)
    for j in range(m):
        prob += (pulp.lpSum([H[j, i] * c[i] for i in range(N)]) == 2 * y[j])
    prob += (pulp.lpSum([dual_logical[i] * c[i] for i in range(N)]) == 1 + 2 * w)
    
    prob.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[prob.status] == 'Optimal':
        return int(pulp.value(prob.objective))
    return float('inf')

def main():
    code_id = "D3_6"
    N = 168
    print(f"--- Running MILP Distance Calculation for {code_id} (N={N}) ---")
    start_time = time.time()

    HX = parse_matrix(f"{code_id}_HX.txt", N)
    HZ = parse_matrix(f"{code_id}_HZ.txt", N)

    log_X = get_logical_basis(HX, HZ, N)
    log_Z = get_logical_basis(HZ, HX, N)

    print(f"Extracted {len(log_X)} logical X and {len(log_Z)} logical Z basis operators.")

    print("Solving MILP for Z-distance (using X-logical operators)...")
    d_Z = min([solve_milp(HX, lx, N) for lx in log_X], default=float('inf'))
    print(f"  -> d_Z = {d_Z}")

    print("Solving MILP for X-distance (using Z-logical operators)...")
    d_X = min([solve_milp(HZ, lz, N) for lz in log_Z], default=float('inf'))
    print(f"  -> d_X = {d_X}")

    exact_d = min(d_Z, d_X)
    elapsed = time.time() - start_time
    print(f"\nFINAL RESULT: {code_id} Exact Distance = {exact_d} (Completed in {elapsed:.2f} seconds)")

if __name__ == "__main__":
    main()

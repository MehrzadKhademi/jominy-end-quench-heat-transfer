import numpy as np
import time

# -----------------------------
# Timer
# -----------------------------
t_start = time.perf_counter()

# -----------------------------
# Parameters (CURRENT WORKFLOW)
# -----------------------------
L = 0.05            # [m]
dx = 1e-5           # [m] 10 microns
dt = 0.01           # [s] 10 ms
t_end = 60.0        # [s]

rho = 7800.0        # [kg/m^3]
cp  = 600.0         # [J/(kg*K)]
k   = 30.0          # [W/(m*K)]
alpha = k / (rho * cp)

h = 12000.0         # [W/(m^2*K)]
T_inf = 25.0        # [°C]

T_init = 740.0      # [°C]

# -----------------------------
# Grid
# -----------------------------
nx = int(L / dx) + 1
x = np.linspace(0.0, L, nx)

nt = int(t_end / dt) + 1
t = np.linspace(0.0, t_end, nt)

Fo = alpha * dt / dx**2
print(f"nx={nx}, nt={nt}, Fo={Fo:.6f}")

# -----------------------------
# Build tridiagonal A for: A * T^{n+1} = RHS
# -----------------------------
a = np.zeros(nx)  # subdiag
b = np.zeros(nx)  # diag
c = np.zeros(nx)  # superdiag

k_over_dx = k / dx

# x=0 boundary (Option B, no ghost):
# (k/dx + h) T0 - (k/dx) T1 = h T_inf
b[0] = k_over_dx + h
c[0] = -k_over_dx
a[0] = 0.0

# interior i=1..nx-2 (Backward Euler + central space)
for i in range(1, nx - 1):
    a[i] = -Fo
    b[i] = 1.0 + 2.0 * Fo
    c[i] = -Fo

# x=L insulated:
# T_N - T_{N-1} = 0  ->  (-1)*T_{N-1} + 1*T_N = 0
a[nx - 1] = -1.0
b[nx - 1] = 1.0
c[nx - 1] = 0.0

# -----------------------------
# Thomas solver
# -----------------------------
def thomas_solve(a, b, c, d):
    n = len(d)
    cp = np.zeros(n)
    dp = np.zeros(n)

    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * cp[i - 1]
        if i < n - 1:
            cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom

    xsol = np.zeros(n)
    xsol[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        xsol[i] = dp[i] - cp[i] * xsol[i + 1]

    return xsol

# -----------------------------
# Time stepping (store full history for animation)
# -----------------------------
T_hist = np.zeros((nt, nx), dtype=np.float32)  # float32 saves disk & RAM
T_curr = np.full(nx, T_init, dtype=np.float32)
T_hist[0, :] = T_curr

RHS = np.zeros(nx, dtype=np.float32)

for n in range(1, nt):
    RHS[0] = h * T_inf
    RHS[1:nx-1] = T_curr[1:nx-1]
    RHS[nx-1] = 0.0

    T_next = thomas_solve(a, b, c, RHS).astype(np.float32)
    T_curr = T_next
    T_hist[n, :] = T_curr

    # optional progress print
    if n % 500 == 0:
        print(f"step {n}/{nt}")

# -----------------------------
# Save new .npy files (overwrite old ones)
# -----------------------------
np.save("x.npy", x)
np.save("t.npy", t)
np.save("T_hist.npy", T_hist)

t_end_run = time.perf_counter()
print(f"Saved x.npy, t.npy, T_hist.npy")
print(f"Runtime: {t_end_run - t_start:.2f} s")

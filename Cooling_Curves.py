# -*- coding: utf-8 -*-
"""
Cooling curves along a 1D Jominy-style bar (simplified heat transfer model)

This script solves transient 1D heat conduction in a bar of length L with:
- Convection boundary at x = 0 (water jet):  -k dT/dx = h (T - T_inf)
- Insulated boundary at x = L:              dT/dx = 0

Numerics:
- Backward Euler in time (unconditionally stable)
- Second-order central difference in space
- Tridiagonal linear system solved by Thomas algorithm

Outputs:
- Cooling curves recorded along the bar (default: every 1 mm)
- A set of highlighted curves at specific distances from the quenched end

Note:
This is intentionally simplified: constant material properties and a constant
"effective" convection coefficient h.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# User-configurable parameters
# =============================================================================

# Geometry
L = 0.05  # [m] bar length (50 mm)

# Mesh (FIXED: use a realistic 1D mesh for a 50 mm bar)
# Old code used dx = 1e-5 m (10 µm), which is computationally excessive here.
dx = 1e-4  # [m] 0.1 mm  -> nx ≈ 501 nodes

# Time
dt = 0.01     # [s]
t_end = 60.0  # [s]

# Material properties (constant)
rho = 7800.0   # [kg/m^3]
cp = 600.0     # [J/(kg*K)]
k = 30.0       # [W/(m*K)]
alpha = k / (rho * cp)

# Boundary conditions
h = 12000.0    # [W/(m^2*K)] effective convection coefficient at x=0
T_inf = 25.0   # [°C] coolant temperature

# Initial condition
T_init = 740.0  # [°C] initial uniform temperature

# Recording density (distance-based, so it remains meaningful if dx changes)
REC_EVERY_MM = 1.0  # record one curve every 1 mm along the bar

# Highlight specific distances from quenched end (mm)
HIGHLIGHT_MM = [0.0, 0.1, 0.5, 1.0, 2.5, 50.0]


# =============================================================================
# Tridiagonal solver (Thomas algorithm)
# =============================================================================
def thomas_solve(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """
    Solve a tridiagonal linear system Ax = d using the Thomas algorithm.

    The system is defined by vectors a, b, c (length n):
      a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    where a[0] is unused and c[n-1] is unused.
    """
    n = len(d)
    cp = np.zeros(n, dtype=float)
    dp = np.zeros(n, dtype=float)

    # Forward sweep
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * cp[i - 1]
        if i < n - 1:
            cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom

    # Back substitution
    x = np.zeros(n, dtype=float)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


# =============================================================================
# Main simulation
# =============================================================================
def main() -> None:
    # Grid
    nx = int(round(L / dx)) + 1
    x = np.linspace(0.0, L, nx)

    nt = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, nt)

    # Dimensionless group (useful for diagnostics; BE is stable for any Fo)
    Fo = alpha * dt / dx**2
    print(f"nx={nx}, nt={nt}, dx={dx:.2e} m, dt={dt:.2e} s, Fo={Fo:.6f}")

    # Record nodes at fixed spacing in mm (robust against dx changes)
    rec_step = max(1, int(round((REC_EVERY_MM * 1e-3) / dx)))
    rec_nodes = np.arange(0, nx, rec_step, dtype=int)
    if rec_nodes[-1] != nx - 1:
        rec_nodes = np.append(rec_nodes, nx - 1)

    T_rec = np.zeros((nt, len(rec_nodes)), dtype=float)

    # Build constant tridiagonal system: A * T^{n+1} = RHS
    a = np.zeros(nx, dtype=float)  # sub-diagonal
    b = np.zeros(nx, dtype=float)  # diagonal
    c = np.zeros(nx, dtype=float)  # super-diagonal

    # Boundary at x=0: convection (enforced at n+1)
    # (k/dx + h) T0 - (k/dx) T1 = h T_inf
    k_over_dx = k / dx
    b[0] = k_over_dx + h
    c[0] = -k_over_dx
    a[0] = 0.0

    # Interior: Backward Euler
    # -Fo*T_{i-1} + (1+2Fo)*T_i - Fo*T_{i+1} = T_i^n
    for i in range(1, nx - 1):
        a[i] = -Fo
        b[i] = 1.0 + 2.0 * Fo
        c[i] = -Fo

    # Boundary at x=L: insulated (Neumann 0), enforced at n+1
    # T_N - T_{N-1} = 0  ->  -1*T_{N-1} + 1*T_N = 0
    a[nx - 1] = -1.0
    b[nx - 1] = 1.0
    c[nx - 1] = 0.0

    # Time stepping
    T_curr = np.full(nx, T_init, dtype=float)
    T_rec[0, :] = T_curr[rec_nodes]

    RHS = np.zeros(nx, dtype=float)
    for n in range(1, nt):
        # x=0 convection RHS
        RHS[0] = h * T_inf

        # interior RHS = previous step temperature
        RHS[1:nx - 1] = T_curr[1:nx - 1]

        # x=L insulated RHS
        RHS[nx - 1] = 0.0

        T_curr = thomas_solve(a, b, c, RHS)
        T_rec[n, :] = T_curr[rec_nodes]

    # Plot
    plt.figure(figsize=(10, 6))

    # log scale cannot show t=0
    t_plot = t[1:]
    for j in range(len(rec_nodes)):
        plt.semilogx(t_plot, T_rec[1:, j], alpha=0.25, color="gray")

    # Highlight representative positions
    for mm in HIGHLIGHT_MM:
        idx = int(round((mm * 1e-3) / dx))
        idx = min(max(idx, 0), nx - 1)
        j = int(np.argmin(np.abs(rec_nodes - idx)))
        plt.semilogx(
            t_plot,
            T_rec[1:, j],
            linewidth=2.0,
            label=f"x = {x[rec_nodes[j]]*1000:.2f} mm",
        )

    plt.xlabel("Time [s] (log scale)")
    plt.ylabel("Temperature [°C]")
    plt.ylim(0, 900)
    plt.xlim(1e-2, 1e2)
    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.title(
        "Cooling curves (1D heat conduction, Backward Euler)\n"
        "Convection at x=0, insulated at x=L"
    )
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()


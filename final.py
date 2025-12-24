import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Utility
# ===============================
def ensure_array(y):
    y = np.array(y)
    if y.ndim == 0:
        y = y.reshape(1)
    return y

# ===============================
# General ODE Solver (Euler + RK4)
# ===============================
class ODESolver:
    """
    General IVP solver supporting Euler and RK4.
    f(t, y): derivative function (scalar or vector)
    """

    def __init__(self, f, t0, y0, dtype=None):
        self.f = f
        self.t0 = t0
        self.y0 = np.array(y0, dtype=dtype)
        if self.y0.ndim == 0:
            self.y0 = self.y0.reshape(1)
        self.dim = self.y0.size

    def _step_euler(self, t, y, dt):
        return y + dt * np.array(self.f(t, y))

    def _step_rk4(self, t, y, dt):
        f = self.f
        k1 = np.array(f(t, y))
        k2 = np.array(f(t + dt/2, y + dt*k1/2))
        k3 = np.array(f(t + dt/2, y + dt*k2/2))
        k4 = np.array(f(t + dt,   y + dt*k3))
        return y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

    def solve(self, t_end=None, dt=None, t_eval=None, method="rk4"):
        if t_eval is not None:
            t_eval = np.asarray(t_eval)
        else:
            t_eval = np.arange(self.t0, t_end+dt, dt)

        y = self.y0.copy()
        ts = [t_eval[0]]
        ys = [y.copy()]

        for i in range(1, len(t_eval)):
            t_prev = t_eval[i-1]
            t_now = t_eval[i]
            dt_now = t_now - t_prev

            if method.lower() == "euler":
                y = self._step_euler(t_prev, y, dt_now)
            else:
                y = self._step_rk4(t_prev, y, dt_now)

            ts.append(t_now)
            ys.append(y.copy())

        return np.array(ts), np.vstack(ys)

# ===============================
# Convert n-th order ODE to 1st-order system
# ===============================
def nth_order_to_system(g, n):
    """
    y^(n) = g(t, y, y', ..., y^(n-1))
    converts to first-order system.
    """
    def F(t, Y):
        Y = np.asarray(Y)
        out = np.zeros_like(Y)
        out[:-1] = Y[1:]      # y' = y1, y1' = y2, ...
        out[-1] = g(t, *Y)    # last item = y^(n)
        return out
    return F

# ===============================
# Example: 2D Projectile motion
# ===============================
g = 9.81
def projectile_system(t, Y):
    x, vx, y, vy = Y
    return np.array([vx, 0, vy, -g])  # dx=vx, dvx=0, dy=vy, dvy=-g

# Initial conditions
speed = 30.0
angle = np.radians(45)
vx0 = speed * np.cos(angle)
vy0 = speed * np.sin(angle)
Y0 = [0, vx0, 0, vy0]

t0 = 0
t_end = 2 * vy0 / g * 1.05
dt = 0.01

solver_rk4 = ODESolver(projectile_system, t0, Y0)
solver_eu  = ODESolver(projectile_system, t0, Y0)

t_rk4, y_rk4 = solver_rk4.solve(t_end=t_end, dt=dt, method="rk4")
t_eu,  y_eu  = solver_eu.solve(t_end=t_end, dt=dt, method="euler")

# Extract trajectories
x_rk4, y_rk4_pos = y_rk4[:,0], y_rk4[:,2]
x_eu,  y_eu_pos  = y_eu[:,0],  y_eu[:,2]

# Analytical solution
x_true = vx0 * t_rk4
y_true = vy0 * t_rk4 - 0.5 * g * t_rk4**2

# ===============================
# Plot
# ===============================
plt.figure(figsize=(8, 5))
plt.plot(x_true, y_true, label="Analytical")
plt.plot(x_rk4, y_rk4_pos, label="RK4")
plt.plot(x_eu, y_eu_pos, label="Euler")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Projectile Motion: Analytical vs RK4 vs Euler")
plt.legend()
plt.grid(True)
plt.show()

def print_ode_solution(ts, ys, labels=None, precision=3, step=10):
    """
    ts: time array
    ys: solution array (shape: [len(ts), dim])
    labels: list of variable names
    precision: number of decimal digits
    step: print every `step` points
    """
    ts = np.asarray(ts)
    ys = np.asarray(ys)

    dim = ys.shape[1] if ys.ndim > 1 else 1
    if labels is None:
        labels = [f"y{i}" for i in range(dim)]
    
    # Print header
    header = "t".ljust(8) + "".join([lbl.rjust(10) for lbl in labels])
    print(header)
    print("-" * len(header))
    
    # Print rows with step
    for t, y in zip(ts[::step], ys[::step]):
        if dim == 1:
            row = f"{t:{8}.{precision}f}{y:{10}.{precision}f}"
        else:
            row = f"{t:{8}.{precision}f}" + "".join([f"{yi:{10}.{precision}f}" for yi in y])
        print(row)


print("\nRK4 Solution (partial, every 10 points):")
print_ode_solution(t_rk4, y_rk4, labels=["x", "vx", "y", "vy"], step=10)

print("\nEuler Solution (partial, every 10 points):")
print_ode_solution(t_eu, y_eu, labels=["x", "vx", "y", "vy"], step=10)


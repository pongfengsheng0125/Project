import numpy as np

# 建立近似解函數 y = y(x)
def make_solution_function(vx, vy):
    vx = np.asarray(vx)
    vy = np.asarray(vy)

    def y_func(xq):
        xq = np.asarray(xq)

        if vy.ndim == 1:
            return np.interp(xq, vx, vy)

        yq = []
        for i in range(vy.shape[0]):
            yq.append(np.interp(xq, vx, vy[i,:]))
        return np.array(yq)
    
    return y_func

# Euler Method 
def euler(f, a=0, b=5, n=100, ya=0):
    h = (b-a)/n

    if np.isscalar(ya):
        vy = np.zeros(n+1, dtype=complex if np.iscomplexobj(ya) else float)
    else:
        m = len(ya)
        vy = np.zeros((m, n+1), dtype=complex if np.iscomplexobj(ya) else float)

    vx = np.zeros(n+1)
    x = a
    y = np.array(ya, copy=True)
    
    vx[0] = x
    if y.ndim > 0:
        vy[:,0] = y
    else:
        vy[0] = y

    for i in range(n):
        y = y + h * f(x, y)
        x = a + (i+1)*h
        vx[i+1] = x
        if vy.ndim == 1:
            vy[i+1] = y
        else:
            vy[:,i+1] = y

    y_func = make_solution_function(vx, vy)
    return vx, vy, y_func

# Runge-Kutta 4 
def rk4(f, a=0, b=5, n=100, ya=0):
    h = (b-a)/n

    if np.isscalar(ya):
        vy = np.zeros(n+1, dtype=complex if np.iscomplexobj(ya) else float)
    else:
        m = len(ya)
        vy = np.zeros((m, n+1), dtype=complex if np.iscomplexobj(ya) else float)

    vx = np.zeros(n+1)
    x = a
    y = np.array(ya, copy=True)

    vx[0] = x
    if y.ndim > 0:
        vy[:,0] = y
    else:
        vy[0] = y

    for i in range(n):
        if np.isscalar(y):
            k1 = h * f(x, y)
            k2 = h * f(x + h/2, y + k1/2)
            k3 = h * f(x + h/2, y + k2/2)
            k4 = h * f(x + h, y + k3)
            y = y + (k1 + 2*k2 + 2*k3 + k4)/6
        else:
            k1 = h * f(x, y)
            k2 = h * f(x + h/2, y + k1/2)
            k3 = h * f(x + h/2, y + k2/2)
            k4 = h * f(x + h, y + k3)
            y = y + (k1 + 2*k2 + 2*k3 + k4)/6

        x = a + (i+1)*h
        vx[i+1] = x
        if vy.ndim == 1:
            vy[i+1] = y
        else:
            vy[:,i+1] = y

    y_func = make_solution_function(vx, vy)
    return vx, vy, y_func

#範例 1：拋物線運動（重力加速度 g = 9.8 m/s²）
#運動方程式：v' = -g, y' = v
#初值：y(0)=0, v(0)=20 m/s

g = 9.8
def f_parabola(t, Y):
    y, v = Y
    return np.array([v, -g])

vx, vy, Y_func = rk4(f_parabola, 0, 5, 100, np.array([0.0, 20.0]))
print("位置 y(2s) =", Y_func(2)[0])
print("速度 v(2s) =", Y_func(2)[1])




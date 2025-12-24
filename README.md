# Python ODE Numerical Solver
(Euler Method & Runge–Kutta Method)

## 1. 程式的原理與功能

常微分方程 (ODE) 的數值解法，解法包括 Euler Method 與 Runge–Kutta Method (RK4)，可處理：

1.單變數 ODE 或 系統 ODE（向量形式）

2.複數或實數初值,

3.解析解未知時的近似解

一般問題形式為：

    dy/dx = f(x, y)
    y(a) = ya

其中：
- x 是自變數
- y 是未知函數
- f(x, y) 為已知的變化率函數
- 已知在 x = a 時，y 的初始值為 ya

對系統 ODE：
    
    Y' = F(x, Y)
    Y(a) = Ya   # Y 為向量

---

## 初值問題與數值解概念

對多數微分方程來說，無法寫出解析解，因此改用數值方法：

- 將區間 [a, b] 切成 n 小段
- 每一小段長度為 h = (b - a) / n
- 逐步從 x = a 計算 y 的近似值

---

## Euler Method (歐拉法)

Euler 法的基本想法是：

    下一個 y ≈ 現在的 y + (步長 × 當前斜率)

對應的數學關係為：

    y_{k+1} = y_k + h * f(x_k, y_k)

其中：
- h 是步長
- f(x_k, y_k) 是在目前點的斜率

特點：
- 實作簡單
- 計算速度快
- 精度較低

---

## Runge–Kutta Method (RK4)

RK4 方法透過在一個步長內，多次估計斜率來提高準確度。

計算過程如下：

    k1 = h * f(x_k, y_k)
    k2 = h * f(x_k + h/2, y_k + k1/2)
    k3 = h * f(x_k + h/2, y_k + k2/2)
    k4 = h * f(x_k + h,   y_k + k3)

更新公式為：

    y_{k+1} = y_k + (k1 + 2*k2 + 2*k3 + k4) / 6

特點：
- 四階精度
- 不需要高階導數
- 適合線性與非線性問題

---

## 高階微分方程與系統表示法

高階微分方程可以轉換為一階系統。

例如：

    y'' = -y

可以改寫為：

    y1 = y
    y2 = y'
    y1' = y2
    y2' = -y1

本程式的 `rk4_system()` 即用來處理此類問題。

---

## 解函數 y(x) 的建立方式

數值解只會得到一組離散點：

    (x0, y0), (x1, y1), (x2, y2), ...

為了讓使用者能在任意 x 上估計 y，本專案建立一個近似解函數：

    y = y(x)

其作法為：

- 使用已計算出的數值點
- 進行分段線性插值
- 在相鄰兩點之間用直線估算 y

此函數可以像一般函數一樣呼叫：

    y_value = y(3.5)

---

## 2. 使用方式

基本步驟：

1.定義微分方程

ODE 必須寫成 dy/dx = f(x, y) 的形式

    def f(x, y):
        return -2*y + np.sin(x)

2.呼叫數值解函數

Euler 方法：

    vx, vy, y_func = euler(f, a=0, b=10, n=100, ya=1.0)

RK4 方法：

    vx, vy, y_func = rk4(f, a=0, b=10, n=100, ya=1.0)

3.使用解函數 y(x)

單個 x 值：

    y_value = y_func(3.5)
    print("y(3.5) =", y_value)

多個 x 值：

    x_vals = [1.0, 2.0, 3.0]
    y_vals = y_func(x_vals)
    print("y([1,2,3]) =", y_vals)

4.範例

範例 1：拋物線運動（重力加速度 g = 9.8 m/s²）

運動方程式：y'' = -g
改寫成：v' = -g, y' = v
初值：y(0)=0, v(0)=20 m/s

    g = 9.8
    def f_parabola(t, Y):
        y, v = Y
        return np.array([v, -g])

    vx, vy, Y_func = rk4(f_parabola, 0, 5, 100, np.array([0.0, 20.0]))
    print("位置 y(2s) =", Y_func(2)[0])
    print("速度 v(2s) =", Y_func(2)[1])

範例 2：指數衰減

ODE：dy/dt = -k * y
初值：y(0) = y0, k > 0

    def f_decay(t, y):
        return -0.5*y

    vx, vy, y_func = rk4(f_decay, 0, 10, 100, 10)

    print("y(5) =", y_func(5))
    print("y([1,2,3]) =", y_func([1,2,3]))

範例 3：簡單諧振動
ODE：y'' + w^2 * y = 0
轉為系統：y1 = y, y2 = y'

    omega = 2.0  # rad/s
    def f_shm(t, Y):
        y, v = Y
        return np.array([v, -omega**2 * y])

    Y0 = np.array([1.0, 0.0])   
    vx, vy, Y_func = rk4(f_shm, 0, 10, 1000, Y0)
    print("y(1s) =", Y_func(1.0)[0])
    print("v(1s) =", Y_func(1.0)[1])

範例 4：阻尼振動

ODE：y'' + 2γ y' + ω² y = 0
系統形式：y1=y, y2=y'
初值：y(0)=1, v(0)=0

    gamma = 0.1
    omega = 2.0

    def f_damped(t, Y):
        y, v = Y
        return np.array([v, -2*gamma*v - omega**2*y])

    Y0 = np.array([1.0, 0.0])
    vx, vy, Y_func = rk4(f_damped, 0, 10, 1000, Y0)

    print("y(5s) =", Y_func(5.0)[0])
    print("v(5s) =", Y_func(5.0)[1])

範例 5：複數 ODE

ODE：dz/dt = i * z, z(0) = 1
解為旋轉的複數軌跡

    def f_complex(t, z):
        return 1j * z

    vx, vy, z_func = rk4(f_complex, 0, 10, 1000, 1.0 + 0j)

    print("z(5) =", z_func(5))

---

## 3. 程式架構

    ├── euler()          # Euler method
    ├── rk4()            # Runge–Kutta (Order 4)
    └── make_solution_function()

---

## 4. 開發過程

1.在選擇題目的時候，沒有什麼想法就在老師提供的範例裡面去了解，剛好自己於對微積分與物理感興趣，尤其最近剛接觸微分方程，對這個的概念感到非常新奇，同時也給我造成很大的困擾，所以就想用程式來解決微分方程。

2.在開發過程中，我最初以拋物線運動作為案例進行實作。但早期的程式僅支援單變數 ODE，例如衰減運動或簡單指數增長，無法處理高阶或多變數系統。為了解決這個限制，我研究了如何將高阶微分方程轉換為一階系統，並改進程式以支援系統 ODE，最終可以處理拋物線運動、簡諧運動等案例。

3.最初程式僅輸出離散數值點，沒有提供可在任意 x 查詢的解函數。在開發過程中，我突然想到許多物理問題需要計算特定時間點的位移或速度，因此將解函數 y(x)加入程式，讓使用者可以直接呼叫函數獲得近似解，而不需自己手動插值。

---

## 5. 參考資料來源

https://github.com/cfgnunes/numerical-methods-python.git ode.py

ChatGPT

---

## 6. 程式修改或增強的內容

在開發時，我參考了現有的 Python ODE 範例程式，其中包括 Euler 法、Taylor 法、RK4 法以及多變數系統的 RK4 實作。

1.為了簡化使用方式，移除 Taylor 高階求導的需求，保留 Euler 與 RK4 方法即可應對大部分初值問題。

2.支援複數運算，在原有程式上加入自動偵測初值是否為複數，並在計算與插值過程中正確處理複數型態。

3.加入解函數 y(x) 支援，透過分段線性插值，使用者可直接呼叫函數獲得近似解。

4.Euler 與 RK4 方法加入多變數支援，可解系統 ODE。
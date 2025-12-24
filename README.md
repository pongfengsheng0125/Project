# Python ODE Numerical Solver  
（Euler Method & Runge–Kutta Method）

## 1. 專案簡介（Project Overview）

本專案實作一套 **Python 常微分方程（ODE）通用數值解器**，用於近似求解初值問題（Initial Value Problem, IVP）：

\[
\frac{dy}{dx} = f(x, y), \quad y(a) = y_a
\]

程式實作多種數值方法，包括：

- Euler Method
- Taylor Method（二階、四階）
- Runge–Kutta Method（RK4）
- Runge–Kutta Method for systems（高階微分方程轉為一階系統）

此外，本專案不僅回傳離散數值解，**進一步建構一個近似的解函數**
\[
y = y(x)
\]
使得使用者可以在任意 \(x\) 上評估數值解。

---

## 2. 數學與數值方法原理（Mathematical Background）

### 2.1 初值問題（IVP）

一般的一階常微分方程初值問題表示為：

\[
\frac{dy}{dx} = f(x, y), \quad y(x_0) = y_0
\]

對於高階微分方程：
\[
y^{(n)} = F(x, y, y', \dots)
\]
可藉由引入變數
\[
y_1 = y,\; y_2 = y',\; \dots
\]
轉換為一階微分方程系統。

---

### 2.2 Euler Method

Euler 法以泰勒展開的一階近似為基礎：

\[
y_{n+1} = y_n + h f(x_n, y_n)
\]

- 方法階數：一階
- 全域誤差：\(O(h)\)
- 優點：概念簡單、實作容易
- 缺點：精度較低、穩定性有限

---

### 2.3 Taylor Method（二階與四階）

利用泰勒展開：

\[
y(x+h) = y(x) + h y'(x) + \frac{h^2}{2} y''(x) + \cdots
\]

本程式實作：
- 二階 Taylor Method
- 四階 Taylor Method（需提供導數函數）

---

### 2.4 Runge–Kutta Method（RK4）

RK4 為常用的高精度顯式方法，其更新公式為：

\[
\begin{aligned}
k_1 &= h f(x_n, y_n) \\
k_2 &= h f(x_n + \tfrac{h}{2}, y_n + \tfrac{k_1}{2}) \\
k_3 &= h f(x_n + \tfrac{h}{2}, y_n + \tfrac{k_2}{2}) \\
k_4 &= h f(x_n + h, y_n + k_3)
\end{aligned}
\]

\[
y_{n+1} = y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\]

- 方法階數：四階
- 全域誤差：\(O(h^4)\)
- 不需計算高階導數
- 適用於線性與非線性系統

---

## 3. 解函數 \( y = y(x) \) 的建構（Solution Function）

數值方法僅會得到離散解點：

\[
(x_0,y_0), (x_1,y_1), \dots, (x_n,y_n)
\]

為了得到一個可連續評估的近似解，本專案利用 **分段線性插值（piecewise linear interpolation）** 定義：

\[
y(x) \approx \text{Interp}(x;\{x_i,y_i\})
\]

程式中透過 `numpy.interp` 實作，並支援：

- 實數 / 複數解
- 單變量 / 多變量（系統）

最終回傳一個可呼叫的函數：

```python
y_func(x_query)

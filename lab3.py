import numpy as np
import matplotlib.pyplot as plt

def basis(t):
    h0 = 2*t**3 - 3*t**2 + 1
    h1 = -2*t**3 + 3*t**2
    h2 = t**3 - 2*t**2 + t
    h3 = t**3 - t**2
    return h0,h1,h2,h3

def compute_tcb_vectors(P, tension, continuity, bias):
    n = len(P)
    T_plus = []
    T_minus = []
    
    for i in range(n):
        if i == 0 or i == n-1:
            T_plus.append(np.array([0.0, 0.0]))
            T_minus.append(np.array([0.0, 0.0]))
        else:
            delta_prev = P[i] - P[i-1]
            delta_next = P[i+1] - P[i]
            
            t = tension[i]
            c = continuity[i]
            b = bias[i]
            
            T_p = (
                (1 - t) * (1 - b) * (1 - c) / 2 * delta_next +
                (1 - t) * (1 + b) * (1 + c) / 2 * delta_prev
            )
            
            T_m = (
                (1 - t) * (1 - b) * (1 + c) / 2 * delta_next +
                (1 - t) * (1 + b) * (1 - c) / 2 * delta_prev
            )
            
            T_plus.append(T_p)
            T_minus.append(T_m)
    
    return T_plus, T_minus

def hermite_curve(P0, P1, T0, T1):
    t = np.linspace(0, 1, 40)
    h0, h1, h2, h3 = basis(t)
    x = h0 * P0[0] + h1 * P1[0] + h2 * T0[0] + h3 * T1[0]
    y = h0 * P0[1] + h1 * P1[1] + h2 * T0[1] + h3 * T1[1]
    return x, y

def generate_tcb_spline(P, tension, continuity, bias):
    T_plus, T_minus = compute_tcb_vectors(P, tension, continuity, bias)
    spline_x = []
    spline_y = []
    
    for i in range(len(P)-1):
        T0 = T_plus[i]
        T1 = T_minus[i+1]
        
        x, y = hermite_curve(P[i], P[i+1], T0, T1)
        spline_x.extend(x)
        spline_y.extend(y)
    
    return spline_x, spline_y

if __name__ == "__main__":
    print("Enter the number of control points (3 or more):")
    n = int(input())
    P=[]
    for i in range(n):
        print(f"Enter coordinates of control point {i+1} (x y):")
        x, y = map(float, input().split())
        P.append(np.array([x, y]))
    
    print("Enter tension, continuity and bias values for each control point:")
    tension = []
    continuity = []
    bias = []
    for i in range(1,n):
        t, c, b = map(float, input(f"Point {i+1} (t c b): ").split())
        tension.append(t)
        continuity.append(c)
        bias.append(b)  
    
    x, y = generate_tcb_spline(P, tension, continuity, bias)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, 'b-', label='TCB-сплайн')
    plt.scatter([p[0] for p in P], [p[1] for p in P], c='red')
    plt.grid(True)
    plt.show()
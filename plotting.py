import matplotlib.pyplot as plt
import numpy as np
import csv

#data = []

u = []
rf = []
with open("C:\\Users\\william\\aba_data\\test.csv", "r") as f:
    data = csv.reader(f)
    next(data)
    for row in data:
        u.append(float(row[0]))
        rf.append(float(row[1]))
#print(data)



RADIUS = 4
def spring_calc_0(l, rho, u, a = 3):
    # l = supporting length
    # d = diameter
    # rho = density
    # u = displacement
    # richard abbot curves
    EC5 = 0.082 * (1 - 0.01 * RADIUS * 2) * rho
    f_h_inter = 0.1253 * rho - 20.32
    k_f_el_0 = 0.1374 * rho - 12.9
    k_f_pl_0 = 0.0047 * rho - 2
    f_h_0 = (k_f_el_0 - k_f_pl_0) * u / ((1 + ((k_f_el_0 - k_f_pl_0) * u / f_h_inter) ** a) ** (1 / a)) + k_f_pl_0 * u
    return f_h_0 * l * RADIUS * 2

print(u)
L = 4
RHO = 650
# Data for plotting

t = np.arange(0.0, 15.0, 0.01)
s = spring_calc_0(L, RHO, t)


fig, ax = plt.subplots()
ax.plot(u, rf)


ax.set(xlabel='Deflection (u)', ylabel='Force (N)',
       title='')
ax.grid()

fig.savefig("test.png")
plt.show()


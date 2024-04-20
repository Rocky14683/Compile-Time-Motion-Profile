import matplotlib.pyplot as plt


sv = open("/Users/rockychen/Desktop/ArbitarySCurveProfile/velocity.txt", "r")
sp = open("/Users/rockychen/Desktop/ArbitarySCurveProfile/Position.txt", "r")
sa = open("/Users/rockychen/Desktop/ArbitarySCurveProfile/acceleration.txt", "r")
# np = open("/Users/rockychen/Desktop/Control theory/competition-code/src/objects/motion_profile/Tposition.txt", "r")
# nv = open("/Users/rockychen/Desktop/Control theory/competition-code/src/objects/motion_profile/Tvelocity.txt", "r")

i = 0
t = 0
Sposition = []
Svelocity = []
Sacceleration = []
Stime = []
# Nposition = []
# Nvelocity = []
# Ntime = []

for t in range(-5, 0):
    Sposition.append(0)
    Svelocity.append(0)
    Sacceleration.append(0)
    Stime.append(t / 0.01)


# for t in range(-10, 1):
#     Nposition.append(0)
#     Nvelocity.append(0)
#     Ntime.append(t / 0.01)

while True:
    svline = sv.readline()
    spline = sp.readline()
    saline = sa.readline()

    if not svline:
       break
    Svelocity.append(float(svline))
    Sposition.append(float(spline))
    Sacceleration.append(float(saline))

    Stime.append(i)

    i += 0.01


# while True:
#     npline = np.readline()
#     nvline = nv.readline()

#     if not npline:
#         break

#     Nposition.append(float(npline))
#     Nvelocity.append(float(nvline))

#     Ntime.append(t)

#     t += 0.01


sv.close()
sp.close()
sa.close()
# np.close()
# nv.close()

figure, axis = plt.subplots(2, 2)


axis[0, 0].plot(Stime, Sposition)
axis[0, 0].set_title("S-Curve Position vs. Time")
axis[0, 0].axis([0, 6, -5, 60])

axis[0, 1].plot(Stime, Svelocity)
axis[0, 1].set_title("S-Curve Velocity vs. Time")
axis[0, 1].axis([0, 6, -12, 25])

axis[1, 0].plot(Stime, Sacceleration)
axis[1, 0].set_title("S-Curve Acceleration vs. Time")
axis[1, 0].axis([0, 6, -35, 35])

# axis[1, 0].plot(Ntime, Nposition)
# axis[1, 0].set_title("T-Curve Position vs. Time")
# axis[1, 0].axis([-1, 6, -1, 20])

# axis[1, 1].plot(Ntime, Nvelocity)
# axis[1, 1].set_title("T-Curve Velocity vs. Time")
# axis[1, 1].axis([-1, 6, -1, 10])

plt.show()

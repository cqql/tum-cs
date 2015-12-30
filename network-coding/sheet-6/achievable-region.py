import matplotlib.pyplot as pl

x = [0, 0, 2, 3, 3, 0]
y = [0, 3, 3, 2, 0, 0]
pl.plot(x, y)
pl.fill(x, y, "g", alpha=0.3)

pl.xlabel("r_1")
pl.ylabel("r_2")
pl.axis([-.5, 3.5, -0.5, 3.5])
pl.savefig("achievable-region.png")

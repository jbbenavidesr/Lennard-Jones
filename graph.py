import matplotlib.pyplot as plt

datafile = "datos.dat"

t = []
x = []
with open(datafile) as file:
    for line in file.readlines():
        sep = line.strip().split(" ")

        t.append(float(sep[0]))
        x.append(float(sep[1]))

t_prom = 0
for i in x:
    t_prom += i

t_prom /= len(x)


def norm(x):
    return x/t_prom


x = list(map(norm, x))


plt.plot(t, x)
plt.show()

from numpy import *


a = array(list(range(8)))
b = array(list(range(5)))
c, d = meshgrid(a, b)
print(c.flatten())
print(d.flatten())
e, f = meshgrid(b, a)
print(" ")
print(e.flatten())
print(f.flatten())
print(tile(a, 2))

import matplotlib.pyplot as plt
import math

ndim = 1000
n = ndim*ndim

times= [1, 2, 3, 4]


def gflops(n, time):
	res = math.pow(2*n, 3) / (time * math.pow(10, 9))  
	return res




for t in times:
	t = gflops(n, t)

plt.plot([4, 9, 16, 25], times)
#plt.axis([0, 30, 0, 500])
plt.ylabel('Gflops/s')
plt.xlabel('Threads')
plt.show()

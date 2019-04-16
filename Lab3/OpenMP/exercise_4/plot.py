import matplotlib.pyplot as plt

plt.figure(1)

plt.subplot(221)
plt.plot([1, 2, 4, 8, 16, 32],[0.000110, 0.000196, 0.000286, 0.000468, 0.000887, 0.001781])
plt.ylabel("Time (s)")
plt.xlabel("Threads")
plt.title("N=10")

plt.subplot(222)
plt.plot([1, 2, 4, 8, 16, 32],[ 0.007253,  0.003845, 0.002346, 0.001785, 0.001615, 0.002515])
plt.ylabel("Time (s)")
plt.xlabel("Threads")
plt.title("N=100")

plt.subplot(223)
plt.plot([1, 2, 4, 8, 16, 32],[0.731963, 0.365002, 0.202308, 0.116530, 0.059643, 0.033757])
plt.ylabel("Time (s)")
plt.xlabel("Threads")
plt.title("N=1000")

plt.subplot(224)
plt.plot([1, 2, 4, 8, 16, 32],[72.777791,  36.605300, 20.060659, 11.562295, 5.783481, 2.903739])
plt.ylabel("Time (s)")
plt.xlabel("Threads")
plt.title("N=10000")


plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.35, wspace=0.35)
plt.show()

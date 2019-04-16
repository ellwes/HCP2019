import matplotlib.pyplot as plt

plt.plot([1,2,4,8,16,32], [11.371048,5.816980,3.181078,1.803991,0.938662,0.471003])
plt.ylabel('Time (s)')
plt.xlabel('Threads')
plt.title('Improvement with multiple threads')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#datareal = np.loadtxt('probmatrix2.txt')
#plt.plot(t, datareal)
datamatrix = np.loadtxt('matrix2.txt')
fig, ax = plt.subplots()

def update(i):
    im_normed = datamatrix[200*i:(200*i + 200),:]
    ax.imshow(im_normed)
    ax.set_axis_off()

anim = FuncAnimation(fig, update, frames=320, interval=1000)
plt.show()
#anim.save('animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)


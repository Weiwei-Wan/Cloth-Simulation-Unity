import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("stable_data/im_energy02.txt", sep=" ")
X_im_02 = df1.iloc[0:, 1].astype(float)
X_im_02 = X_im_02.to_numpy()[0:100]

df2 = pd.read_csv("stable_data/im_energy04.txt", sep=" ")
X_im_04 = df2.iloc[0:, 1].astype(float)
X_im_04 = X_im_04.to_numpy()[0:100]

df3 = pd.read_csv("stable_data/im_energy05.txt", sep=" ")
X_im_05 = df3.iloc[0:, 1].astype(float)
X_im_05 = X_im_05.to_numpy()[0:100]

df4 = pd.read_csv("stable_data/im_energy06.txt", sep=" ")
X_im_06 = df4.iloc[0:, 1].astype(float)
X_im_06 = X_im_06.to_numpy()[0:100]

plt.figure(figsize=(6, 8))

x=np.arange(0,100)
l1=plt.plot(x,X_im_02,'y--',label='dt = 0.02s')
l2=plt.plot(x,X_im_04,'r--',label='dt = 0.04s')
l3=plt.plot(x,X_im_05,'g--',label='dt = 0.05s')
l4=plt.plot(x,X_im_06,'b--',label='dt = 0.06s')
plt.plot(x,X_im_02,'y*-', x,X_im_04,'ro-',x,X_im_05,'g+-',x,X_im_06,'b^-', markersize='2')
plt.xlabel('Frame')
plt.ylabel('Spring Energy')
plt.ylim((1e-5, 1e31))
plt.legend()

plt.yscale('log')
plt.savefig('./im_stable.jpg', dpi=600)
plt.show()
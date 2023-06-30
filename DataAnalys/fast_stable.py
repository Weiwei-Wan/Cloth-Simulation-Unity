import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("stable_data/fast_energy002.txt", sep=" ")
X_fast_002 = df1.iloc[0:, 1].astype(float)
X_fast_002 = X_fast_002.to_numpy()[0:100]

df2 = pd.read_csv("stable_data/fast_energy02.txt", sep=" ")
X_fast_02 = df2.iloc[0:, 1].astype(float)
X_fast_02 = X_fast_02.to_numpy()[0:100]

df3 = pd.read_csv("stable_data/fast_energy2.txt", sep=" ")
X_fast_2 = df3.iloc[0:, 1].astype(float)
X_fast_2 = X_fast_2.to_numpy()[0:100]

df4 = pd.read_csv("stable_data/fast_energy1.txt", sep=" ")
X_fast_1 = df4.iloc[0:, 1].astype(float)
X_fast_1 = X_fast_1.to_numpy()[0:100]

plt.figure(figsize=(6, 8))

x=np.arange(0,100)
l1=plt.plot(x,X_fast_002,'y--',label='dt = 0.002s')
l2=plt.plot(x,X_fast_02,'r--',label='dt = 0.02s')
l3=plt.plot(x,X_fast_2,'g--',label='dt = 0.2s')
l4=plt.plot(x,X_fast_1,'b--',label='dt = 1s')
plt.plot(x,X_fast_002,'y*-', x,X_fast_02,'ro-',x,X_fast_2,'g+-',x,X_fast_1,'b^-', markersize='2')
plt.xlabel('Frame')
plt.ylabel('Spring Energy')
plt.ylim((1e-5, 1e31))
plt.legend()

plt.yscale('log')
plt.savefig('./fast_stable.jpg', dpi=600)
plt.show()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("stable_data/ex_energy_002.txt", sep=" ")
X_ex_002 = df1.iloc[0:, 1].astype(float)
X_ex_002 = X_ex_002.to_numpy()[0:100]

df2 = pd.read_csv("stable_data/ex_energy006.txt", sep=" ")
X_ex_006 = df2.iloc[0:, 1].astype(float)
X_ex_006 = X_ex_006.to_numpy()[0:100]

df3 = pd.read_csv("stable_data/ex_energy007.txt", sep=" ")
X_ex_007 = df3.iloc[0:, 1].astype(float)
X_ex_007 = X_ex_007.to_numpy()[0:100]

df4 = pd.read_csv("stable_data/ex_energy008.txt", sep=" ")
X_ex_008 = df4.iloc[0:, 1].astype(float)
X_ex_008 = X_ex_008.to_numpy()[0:100]

plt.figure(figsize=(3.2, 4))

x=np.arange(0,100)
l1=plt.plot(x,X_ex_002,'y--',label='dt = 0.002s')
l2=plt.plot(x,X_ex_006,'r--',label='dt = 0.006s')
l3=plt.plot(x,X_ex_007,'g--',label='dt = 0.007s')
l4=plt.plot(x,X_ex_008,'b--',label='dt = 0.008s')
plt.plot(x,X_ex_002,'y*-', x,X_ex_006,'ro-',x,X_ex_007,'g+-',x,X_ex_008,'b^-', markersize='2')
plt.xlabel('Frame')
plt.ylabel('Spring Energy')
plt.ylim((1e-5, 1e31))
plt.legend()

plt.yscale('log')
plt.savefig('./ex_stable2.jpg', dpi=600)
plt.show()
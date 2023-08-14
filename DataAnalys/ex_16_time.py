import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("ex/ex-16-001.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_001 = df1.iloc[1:, 2]
X_ex_16_001 = df1.iloc[1:, 2].astype(float)
X_ex_16_001 = X_ex_16_001.to_numpy()
X_ex_16_001_Mean = np.mean(X_ex_16_001)
X_ex_16_001_StdErr = np.std(X_ex_16_001, ddof=1) / np.sqrt(np.size(X_ex_16_001))

df2 = pd.read_csv("ex/ex-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_002 = df2.iloc[1:, 2]
X_ex_16_002 = df2.iloc[1:, 2].astype(float)
X_ex_16_002 = X_ex_16_002.to_numpy()
X_ex_16_002_Mean = np.mean(X_ex_16_002)
X_ex_16_002_StdErr = np.std(X_ex_16_002, ddof=1) / np.sqrt(np.size(X_ex_16_002))

df3 = pd.read_csv("ex/ex-16-005.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_005 = df3.iloc[1:, 2]
X_ex_16_005 = df3.iloc[1:, 2].astype(float)
X_ex_16_005 = X_ex_16_005.to_numpy()
X_ex_16_005_Mean = np.mean(X_ex_16_005)
X_ex_16_005_StdErr = np.std(X_ex_16_005, ddof=1) / np.sqrt(np.size(X_ex_16_005))

df4 = pd.read_csv("ex/ex-16-0005.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_0005 = df4.iloc[1:, 2]
X_ex_16_0005 = df4.iloc[1:, 2].astype(float)
X_ex_16_0005 = X_ex_16_0005.to_numpy()
X_ex_16_0005_Mean = np.mean(X_ex_16_0005)
X_ex_16_0005_StdErr = np.std(X_ex_16_0005, ddof=1) / np.sqrt(np.size(X_ex_16_0005))

df5 = pd.read_csv("ex/ex-16-0001.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_0001 = df5.iloc[1:, 2]
X_ex_16_0001 = df5.iloc[1:, 2].astype(float)
X_ex_16_0001 = X_ex_16_0001.to_numpy()
X_ex_16_0001_Mean = np.mean(X_ex_16_0001)
X_ex_16_0001_StdErr = np.std(X_ex_16_0001, ddof=1) / np.sqrt(np.size(X_ex_16_0001))

plt.figure(figsize=(3, 4))

x=[0.0001, 0.0005, 0.001,0.002,0.005]
y=[X_ex_16_0001_Mean, X_ex_16_0005_Mean, X_ex_16_001_Mean, X_ex_16_002_Mean, X_ex_16_005_Mean]
std_err=[X_ex_16_0001_StdErr, X_ex_16_0005_StdErr, X_ex_16_001_StdErr, X_ex_16_002_StdErr, X_ex_16_005_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.errorbar(x,y,fmt='ro-',yerr=std_err)
plt.xlabel('Time step (ms)')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./fig_X_ex_16_time.jpg', dpi=600)
plt.show()
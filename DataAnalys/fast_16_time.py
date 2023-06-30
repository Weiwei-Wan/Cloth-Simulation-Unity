import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("fast/fast-16-001.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_001 = df1.iloc[1:, 2]
X_fast_16_001 = df1.iloc[1:, 2].astype(float)
X_fast_16_001 = X_fast_16_001.to_numpy()
X_fast_16_001_Mean = np.mean(X_fast_16_001)
X_fast_16_001_StdErr = np.std(X_fast_16_001, ddof=1) / np.sqrt(np.size(X_fast_16_001))

df2 = pd.read_csv("fast/fast-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_002 = df2.iloc[1:, 2]
X_fast_16_002 = df2.iloc[1:, 2].astype(float)
X_fast_16_002 = X_fast_16_002.to_numpy()
X_fast_16_002_Mean = np.mean(X_fast_16_002)
X_fast_16_002_StdErr = np.std(X_fast_16_002, ddof=1) / np.sqrt(np.size(X_fast_16_002))

df3 = pd.read_csv("fast/fast-16-005.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_005 = df3.iloc[1:, 2]
X_fast_16_005 = df3.iloc[1:, 2].astype(float)
X_fast_16_005 = X_fast_16_005.to_numpy()
X_fast_16_005_Mean = np.mean(X_fast_16_005)
X_fast_16_005_StdErr = np.std(X_fast_16_005, ddof=1) / np.sqrt(np.size(X_fast_16_005))

df4 = pd.read_csv("fast/fast-16-01.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_01 = df4.iloc[1:, 2]
X_fast_16_01 = df4.iloc[1:, 2].astype(float)
X_fast_16_01 = X_fast_16_01.to_numpy()
X_fast_16_01_Mean = np.mean(X_fast_16_01)
X_fast_16_01_StdErr = np.std(X_fast_16_01, ddof=1) / np.sqrt(np.size(X_fast_16_01))

df5 = pd.read_csv("fast/fast-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_02 = df5.iloc[1:, 2]
X_fast_16_02 = df5.iloc[1:, 2].astype(float)
X_fast_16_02 = X_fast_16_02.to_numpy()
X_fast_16_02_Mean = np.mean(X_fast_16_02)
X_fast_16_02_StdErr = np.std(X_fast_16_02, ddof=1) / np.sqrt(np.size(X_fast_16_02))

plt.figure(figsize=(6, 8))

x=[0.001,0.002,0.005,0.01,0.02]
y=[X_fast_16_001_Mean, X_fast_16_002_Mean, X_fast_16_005_Mean, X_fast_16_01_Mean, X_fast_16_02_Mean]
std_err=[X_fast_16_001_StdErr, X_fast_16_002_StdErr, X_fast_16_005_StdErr, X_fast_16_01_StdErr, X_fast_16_02_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.errorbar(x,y,fmt='ro-',yerr=std_err)
plt.xlabel('Time step (ms)')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./fig_X_fast_16_time.jpg', dpi=600)
plt.show()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("ex/ex-8-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_8_002 = df1.iloc[1:, 2].astype(float)
X_ex_8_002 = X_ex_8_002.to_numpy()
X_ex_8_002_Mean = np.mean(X_ex_8_002)
X_ex_8_002_StdErr = np.std(X_ex_8_002, ddof=1) / np.sqrt(np.size(X_ex_8_002))

df2 = pd.read_csv("ex/ex-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_002 = df2.iloc[1:, 2]
X_ex_16_002 = df2.iloc[1:, 2].astype(float)
X_ex_16_002 = X_ex_16_002.to_numpy()
X_ex_16_002_Mean = np.mean(X_ex_16_002)
X_ex_16_002_StdErr = np.std(X_ex_16_002, ddof=1) / np.sqrt(np.size(X_ex_16_002))

df3 = pd.read_csv("ex/ex-32-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_32_002 = df3.iloc[1:, 2]
X_ex_32_002 = df3.iloc[1:, 2].astype(float)
X_ex_32_002 = X_ex_32_002.to_numpy()
X_ex_32_002_Mean = np.mean(X_ex_32_002)
X_ex_32_002_StdErr = np.std(X_ex_32_002, ddof=1) / np.sqrt(np.size(X_ex_32_002))

df4 = pd.read_csv("ex/ex-64-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_64_002 = df4.iloc[1:, 2]
X_ex_64_002 = df4.iloc[1:, 2].astype(float)
X_ex_64_002 = X_ex_64_002.to_numpy()
X_ex_64_002_Mean = np.mean(X_ex_64_002)
X_ex_64_002_StdErr = np.std(X_ex_64_002, ddof=1) / np.sqrt(np.size(X_ex_64_002))

plt.figure(figsize=(5, 4))

x=[1,2,3]
y=[X_ex_8_002_Mean, X_ex_16_002_Mean, X_ex_32_002_Mean]
std_err=[X_ex_8_002_StdErr, X_ex_16_002_StdErr, X_ex_32_002_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.bar(x,y,color=['b','g','r'],yerr=std_err,error_kw=error_params,tick_label=['8x8','16x16','32x32'])
plt.xlabel('Cloth size')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./fig_X_ex_002_size.jpg', dpi=600)
plt.show()
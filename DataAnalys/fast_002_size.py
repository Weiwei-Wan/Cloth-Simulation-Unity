
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("fast/fast-8-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_8_02 = df1.iloc[1:, 2].astype(float)
X_fast_8_02 = X_fast_8_02.to_numpy()
X_fast_8_02_Mean = np.mean(X_fast_8_02)
X_fast_8_02_StdErr = np.std(X_fast_8_02, ddof=1) / np.sqrt(np.size(X_fast_8_02))

df2 = pd.read_csv("fast/fast-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_02 = df2.iloc[1:, 2]
X_fast_16_02 = df2.iloc[1:, 2].astype(float)
X_fast_16_02 = X_fast_16_02.to_numpy()
X_fast_16_02_Mean = np.mean(X_fast_16_02)
X_fast_16_02_StdErr = np.std(X_fast_16_02, ddof=1) / np.sqrt(np.size(X_fast_16_02))

df3 = pd.read_csv("fast/fast-32-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_32_02 = df3.iloc[1:, 2]
X_fast_32_02 = df3.iloc[1:, 2].astype(float)
X_fast_32_02 = X_fast_32_02.to_numpy()
X_fast_32_02_Mean = np.mean(X_fast_32_02)
X_fast_32_02_StdErr = np.std(X_fast_32_02, ddof=1) / np.sqrt(np.size(X_fast_32_02))

plt.figure(figsize=(3, 4))

x=[1,2,3]
y=[X_fast_8_02_Mean, X_fast_16_02_Mean, X_fast_32_02_Mean]
std_err=[X_fast_8_02_StdErr, X_fast_16_02_StdErr, X_fast_32_02_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.bar(x,y,color=['b','g','r'],yerr=std_err,error_kw=error_params,tick_label=['8x8','16x16','32x32'])
plt.xlabel('Cloth size')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./fig_X_fast_02_size.jpg', dpi=600)
plt.show()
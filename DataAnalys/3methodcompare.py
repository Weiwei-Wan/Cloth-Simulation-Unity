import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("ex/ex-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_16_002 = df1.iloc[1:, 2].astype(float)
X_ex_16_002 = X_ex_16_002.to_numpy()
X_ex_16_002_Mean = np.mean(X_ex_16_002)
X_ex_16_002_StdErr = np.std(X_ex_16_002, ddof=1) / np.sqrt(np.size(X_ex_16_002))

df2 = pd.read_csv("fast/fast-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_16_02 = df2.iloc[1:, 2].astype(float)
X_fast_16_02 = X_fast_16_02.to_numpy()
X_fast_16_02_Mean = np.mean(X_fast_16_02)
X_fast_16_02_StdErr = np.std(X_fast_16_02, ddof=1) / np.sqrt(np.size(X_fast_16_02))

df3 = pd.read_csv("im/im-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_02 = df3.iloc[1:, 2].astype(float)
X_im_16_02 = X_im_16_02.to_numpy()
X_im_16_02_Mean = np.mean(X_im_16_02)
X_im_16_02_StdErr = np.std(X_im_16_02, ddof=1) / np.sqrt(np.size(X_im_16_02))

plt.figure(figsize=(6, 8))

x=[1,2,3]
y=[X_ex_16_002_Mean, X_im_16_02_Mean, X_fast_16_02_Mean]
std_err=[X_ex_16_002_StdErr, X_im_16_02_StdErr, X_fast_16_02_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.bar(x,y,color=['b','g','r'],yerr=std_err,error_kw=error_params,tick_label=['Explicit','Implict','Fast'])
plt.xlabel('Cloth size')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./3_method_compare.jpg', dpi=600)
plt.show()
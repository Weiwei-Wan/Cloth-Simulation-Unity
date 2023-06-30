import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df1 = pd.read_csv("ex/ex-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_cpu = df1.iloc[1:, 2].astype(float)
X_ex_cpu = X_ex_cpu.to_numpy()
X_ex_cpu_Mean = np.mean(X_ex_cpu)
X_ex_cpu_StdErr = np.std(X_ex_cpu, ddof=1) / np.sqrt(np.size(X_ex_cpu))

df2 = pd.read_csv("GPU/ex-gpu-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_ex_gpu = df2.iloc[1:, 2].astype(float)
X_ex_gpu = X_ex_gpu.to_numpy()
X_ex_gpu_Mean = np.mean(X_ex_gpu)
X_ex_gpu_StdErr = np.std(X_ex_gpu, ddof=1) / np.sqrt(np.size(X_ex_gpu))

df3 = pd.read_csv("fast/fast-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_cpu = df3.iloc[1:, 2].astype(float)
X_fast_cpu = X_fast_cpu.to_numpy()
X_fast_cpu_Mean = np.mean(X_fast_cpu)
X_fast_cpu_StdErr = np.std(X_fast_cpu, ddof=1) / np.sqrt(np.size(X_fast_cpu))

df4 = pd.read_csv("GPU/fast-gpu-0.02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_fast_gpu = df4.iloc[1:, 2].astype(float)
X_fast_gpu = X_fast_gpu.to_numpy()
X_fast_gpu_Mean = np.mean(X_fast_gpu)
X_fast_gpu_StdErr = np.std(X_fast_gpu, ddof=1) / np.sqrt(np.size(X_fast_gpu))

df5 = pd.read_csv("im/im-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_cpu = df5.iloc[1:, 2].astype(float)
X_im_cpu = X_im_cpu.to_numpy()
X_im_cpu_Mean = np.mean(X_im_cpu)
X_im_cpu_StdErr = np.std(X_im_cpu, ddof=1) / np.sqrt(np.size(X_im_cpu))

df6 = pd.read_csv("GPU/im-gpu-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_gpu = df6.iloc[1:, 2].astype(float)
X_im_gpu = X_im_gpu.to_numpy()
X_im_gpu_Mean = np.mean(X_im_gpu)
X_im_gpu_StdErr = np.std(X_im_gpu, ddof=1) / np.sqrt(np.size(X_im_gpu))

plt.figure(figsize=(6, 8))

x=[1,2,3]
y1=[X_ex_cpu_Mean, X_im_cpu_Mean, X_fast_cpu_Mean]
y2=[X_ex_gpu_Mean, X_im_gpu_Mean, X_fast_gpu_Mean]
std_err1=[X_ex_cpu_StdErr, X_im_cpu_StdErr, X_fast_cpu_StdErr]
std_err2=[X_ex_gpu_StdErr, X_im_gpu_StdErr, X_fast_gpu_StdErr]

error_params=dict(elinewidth=3,ecolor='k',capsize=5)
a = plt.bar(x,y1,width=0.3,yerr=std_err1,error_kw=error_params, label='CPU', fc="b", tick_label=['Explicit','Implict','Fast'])
x=[1.4,2.4,3.4]
b = plt.bar(x,y2,width=0.3,yerr=std_err2,error_kw=error_params, label='GPU', fc="g")

plt.xlabel('Algorithm')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.legend()
plt.savefig('./cpu_gpu_compare.jpg', dpi=600)
plt.show()
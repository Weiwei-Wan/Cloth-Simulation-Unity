import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read the data
df2 = pd.read_csv("im/im-16-002.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_002 = df2.iloc[1:, 2]
X_im_16_002 = df2.iloc[1:, 2].astype(float)
X_im_16_002 = X_im_16_002.to_numpy()
X_im_16_002_Mean = np.mean(X_im_16_002)
X_im_16_002_StdErr = np.std(X_im_16_002, ddof=1) / np.sqrt(np.size(X_im_16_002))

df1 = pd.read_csv("im/im-16-004.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_004 = df1.iloc[1:, 2]
X_im_16_004 = df1.iloc[1:, 2].astype(float)
X_im_16_004 = X_im_16_004.to_numpy()
X_im_16_004_Mean = np.mean(X_im_16_004)
X_im_16_004_StdErr = np.std(X_im_16_004, ddof=1) / np.sqrt(np.size(X_im_16_004))

df3 = pd.read_csv("im/im-16-006.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_006 = df3.iloc[1:, 2]
X_im_16_006 = df3.iloc[1:, 2].astype(float)
X_im_16_006 = X_im_16_006.to_numpy()
X_im_16_006_Mean = np.mean(X_im_16_006)
X_im_16_006_StdErr = np.std(X_im_16_006, ddof=1) / np.sqrt(np.size(X_im_16_006))

df4 = pd.read_csv("im/im-16-008.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_008 = df4.iloc[1:, 2]
X_im_16_008 = df4.iloc[1:, 2].astype(float)
X_im_16_008 = X_im_16_008.to_numpy()
X_im_16_008_Mean = np.mean(X_im_16_008)
X_im_16_008_StdErr = np.std(X_im_16_008, ddof=1) / np.sqrt(np.size(X_im_16_008))

df5 = pd.read_csv("im/im-16-01.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_01 = df5.iloc[1:, 2]
X_im_16_01 = df5.iloc[1:, 2].astype(float)
X_im_16_01 = X_im_16_01.to_numpy()
X_im_16_01_Mean = np.mean(X_im_16_01)
X_im_16_01_StdErr = np.std(X_im_16_01, ddof=1) / np.sqrt(np.size(X_im_16_01))

df6 = pd.read_csv("im/im-16-012.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_012 = df6.iloc[1:, 2]
X_im_16_012 = df6.iloc[1:, 2].astype(float)
X_im_16_012 = X_im_16_012.to_numpy()
X_im_16_012_Mean = np.mean(X_im_16_012)
X_im_16_012_StdErr = np.std(X_im_16_012, ddof=1) / np.sqrt(np.size(X_im_16_012))

df7 = pd.read_csv("im/im-16-014.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_014 = df7.iloc[1:, 2]
X_im_16_014 = df7.iloc[1:, 2].astype(float)
X_im_16_014 = X_im_16_014.to_numpy()
X_im_16_014_Mean = np.mean(X_im_16_014)
X_im_16_014_StdErr = np.std(X_im_16_014, ddof=1) / np.sqrt(np.size(X_im_16_014))

df8 = pd.read_csv("im/im-16-016.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_016 = df8.iloc[1:, 2]
X_im_16_016 = df8.iloc[1:, 2].astype(float)
X_im_16_016 = X_im_16_016.to_numpy()
X_im_16_016_Mean = np.mean(X_im_16_016)
X_im_16_016_StdErr = np.std(X_im_16_016, ddof=1) / np.sqrt(np.size(X_im_16_016))

df9 = pd.read_csv("im/im-16-018.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_018 = df9.iloc[1:, 2]
X_im_16_018 = df9.iloc[1:, 2].astype(float)
X_im_16_018 = X_im_16_018.to_numpy()
X_im_16_018_Mean = np.mean(X_im_16_018)
X_im_16_018_StdErr = np.std(X_im_16_018, ddof=1) / np.sqrt(np.size(X_im_16_018))

df10 = pd.read_csv("im/im-16-02.csv", sep=",", names=["a", "b", "frame_time", "c"])
X_im_16_02 = df10.iloc[1:, 2]
X_im_16_02 = df10.iloc[1:, 2].astype(float)
X_im_16_02 = X_im_16_02.to_numpy()
X_im_16_02_Mean = np.mean(X_im_16_02)
X_im_16_02_StdErr = np.std(X_im_16_02, ddof=1) / np.sqrt(np.size(X_im_16_02))

plt.figure(figsize=(3, 4))

x=[0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02]
y=[X_im_16_002_Mean, X_im_16_004_Mean, X_im_16_006_Mean, X_im_16_008_Mean, X_im_16_01_Mean, X_im_16_012_Mean, X_im_16_014_Mean, X_im_16_016_Mean, X_im_16_018_Mean, X_im_16_02_Mean]
std_err=[X_im_16_002_StdErr, X_im_16_004_StdErr, X_im_16_006_StdErr, X_im_16_008_StdErr, X_im_16_01_StdErr, X_im_16_012_StdErr, X_im_16_014_StdErr, X_im_16_016_StdErr, X_im_16_018_StdErr, X_im_16_02_StdErr]
error_params=dict(elinewidth=3,ecolor='k',capsize=5)
plt.errorbar(x,y,fmt='ro-',yerr=std_err)
plt.xlabel('Time step (ms)')
plt.ylabel('Time / Frame (ms)')
plt.yscale('log')
plt.savefig('./fig_X_im_16_time.jpg', dpi=600)
plt.show()
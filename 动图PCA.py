import sys
print(sys.executable)
import pandas as pd
import plotly.express as px
import plotly.io as pio
import numpy as np
from scipy.stats import chi2
import plotly.graph_objects as go

# 读取eigenvec文件
eigenvec_path = r"Z:\临时下载文件\大好河山芯片数据\大好河山芯片数据\pca_output_replaced.eigenvec"
df = pd.read_csv(eigenvec_path, delim_whitespace=True, header=None)

# 假设文件中有 "FID", "IID", "PC1", "PC2", "PC3" 这些列
df.columns = ['FID', 'IID', 'PC1', 'PC2', 'PC3']

print(df.head())
print("数据行数：", len(df))

# 使用Plotly绘制3D散点图
fig = px.scatter_3d(df, x='PC1', y='PC2', z='PC3', color='PC1', 
                    hover_data=['FID', 'IID'], 
                    color_continuous_scale='Viridis', 
                    title="3D PCA Plot")

# 设置背景颜色和图形风格
fig.update_layout(scene=dict(
                    xaxis_title='PC1',
                    yaxis_title='PC2',
                    zaxis_title='PC3',
                    bgcolor='rgba(240,240,240,0.9)'  # 背景颜色
                ))

# 计算均值和协方差
mean = df[['PC1', 'PC2', 'PC3']].mean().values
cov = np.cov(df[['PC1', 'PC2', 'PC3']].values, rowvar=False)

# 置信水平（95%置信区间）
confidence_level = 0.95
radius = np.sqrt(chi2.ppf(confidence_level, 3))

# 生成椭球点云
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones_like(u), np.cos(v))
ellipsoid = np.stack((x, y, z), axis=-1)

# 奇异值分解，得到椭球的轴
U, s, rotation = np.linalg.svd(cov)
axes = radius * np.sqrt(s)

# 变换椭球点云
for i in range(len(u)):
    for j in range(len(v)):
        ellipsoid[i, j, :] = mean + np.dot(ellipsoid[i, j, :] * axes, rotation)

# 添加到plotly图中
ellipsoid_trace = go.Surface(
    x=ellipsoid[:, :, 0],
    y=ellipsoid[:, :, 1],
    z=ellipsoid[:, :, 2],
    opacity=0.2,
    colorscale='Blues',
    showscale=False,
    name='95%置信椭球'
)

fig.add_trace(ellipsoid_trace)

# 保存为HTML文件
pio.write_html(fig, file="pca_plot_with_ellipsoid.html", auto_open=True)

# 如果希望保存为静态图片
# fig.write_image("pca_plot.png")

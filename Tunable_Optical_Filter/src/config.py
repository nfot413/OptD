# -*- coding: utf-8 -*-
"""
光学滤波器优化项目的配置文件

此文件使用简单的全局变量来定义所有参数，便于直接导入和使用。
变量名使用大写，以遵循Python中对常量的通用约定 (PEP 8)。
"""
import numpy as np

# --- 1. 结构参数 ---
# 定义滤波器的物理构成
N_KU = 2  # 上臂的微环数量
M_KL = 1  # 下臂的微环数量

# --- 2. 物理参数 ---
# 定义器件和系统的物理属性
F_CENTER = 193.1e12          # 中心频率 (Hz)
FSR = 100e9                  # 自由光谱范围 (Hz)
T_SELF_COUPLING = 0.979888   # 微环的自耦合/往返传输系数
THETA_I = np.pi / 2          # 输入MZI的内部相移
THETA_O = np.pi / 2          # 输出MZI的内部相移
PHI_T = 0.0                  # 上臂的附加相移
PHI_B = 0.0                  # 下臂的附加相移
DELAY_VAL = 1.0              # 下臂延迟线的时间延迟
DELAY_PHASE = 0.0            # 下臂延迟线的附加相移

# --- 3. 仿真参数 ---
# 控制仿真的范围和精度
S_FSR_SPAN = 10              # 仿真频率范围，中心频率 +/- s*FSR
W_MIN_PI_MULTIPLE = -20      # 归一化角频率范围的下限 (乘以pi)
W_MAX_PI_MULTIPLE = 20       # 归一化角频率范围的上限 (乘以pi)
DW_STEP = 0.006285           # 归一化角频率的步长

# --- 4. 目标滤波器参数 ---
# 定义理想滤波器的形状
BANDWIDTH = 50e9             # 通带宽度 (Hz)
PASSBAND_DB = 0              # 通带电平 (dB)
STOPBAND_DB = -40            # 阻带电平 (dB)

# --- 5. 优化器参数 ---
# 配置差分进化算法
MAX_ITERATIONS = 300         # 最大迭代次数
POPULATION_SIZE = 20         # 种群大小

# --- 6. 绘图参数 ---
# 控制最终输出图像的样式
PLOT_TITLE = "光学滤波器优化结果 (函数式重构)"
PLOT_XLABEL = "频率 (Hz)"
PLOT_YLABEL = "幅度 (dB)"
PLOT_YLIM_MIN = -60
PLOT_YLIM_MAX = 5
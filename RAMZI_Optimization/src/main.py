# src/main.py
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skopt import gp_minimize
from skopt.space import Real
from skopt.plots import plot_convergence

from . import config
from .hardware import OSA, VoltageSource
from .objective import ObjectiveFunction

def run():
    """主优化流程"""
    print("开始进行优化...")

    # --- 初始化硬件 ---
    voltage_source = VoltageSource(config.VOLTAGE_SOURCE_ADDRESS, mock_mode=config.MOCK_HARDWARE)
    osa = OSA(config.OSA_ADDRESS, mock_mode=config.MOCK_HARDWARE)

    # --- 加载目标数据 ---
    try:
        df = pd.read_excel(config.TARGET_SPECTRUM_FILE, header=None)
        target_spectrum = df.iloc[0, 1:1402].values
        print(f"成功加载目标光谱，数据点数: {len(target_spectrum)}")
    except Exception as e:
        print(f"警告: 加载目标文件 '{config.TARGET_SPECTRUM_FILE}' 失败 ({e})。")
        print("将使用内置的默认高斯曲线作为目标。")
        num_points = int((config.OSA_WAVELENGTH_STOP - config.OSA_WAVELENGTH_START) / config.OSA_RESOLUTION) + 1
        target_spectrum = -30 * np.exp(-(np.linspace(-3, 3, num_points))**2)

    # --- 定义优化问题 ---
    dimensions = [Real(name=f'v_{i}', low=config.VOLTAGE_RANGE[0], high=config.VOLTAGE_RANGE[1]) for i in range(config.NUM_VECTORS)]
    objective = ObjectiveFunction(osa, voltage_source, target_spectrum)

    # --- 运行优化 ---
    print(f"启动贝叶斯优化器，总共进行 {config.OPTIMIZER_CALLS} 次测量...")
    start_time = time.time()
    
    result = gp_minimize(
        func=objective,
        dimensions=dimensions,
        n_calls=config.OPTIMIZER_CALLS,
        n_initial_points=config.OPTIMIZER_RANDOM_STARTS,
        random_state=123,
    )
    
    end_time = time.time()
    print(f"\n--- 优化完成 (总耗时: {(end_time - start_time):.2f} 秒) ---")
    print(f"找到的最优电压: {[f'{v:.3f}' for v in result.x]}")
    print(f"达到的最小成本值 (MSE): {result.fun:.6f}")
    
    # --- 绘制结果 ---
    plot_convergence(result)
    plt.show()

if __name__ == '__main__':
    run()
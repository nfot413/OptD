# src/objective.py
import time
import numpy as np
from .hardware import OSA, VoltageSource

class ObjectiveFunction:
    """封装一次完整的硬件评估流程。"""
    def __init__(self, osa: OSA, vs: VoltageSource, target_spectrum: np.ndarray):
        self.osa = osa
        self.voltage_source = vs
        self.target_spectrum = target_spectrum
        self.iteration_count = 0

    def __call__(self, voltages: list[float]) -> float:
        """执行一次评估，返回成本值 (越小越好)。"""
        self.iteration_count += 1
        print(f"\n--- 第 {self.iteration_count} 次评估 ---")
        
        try:
            self.voltage_source.set_all_voltages(voltages)
            time.sleep(0.1)
            _, measured_power = self.osa.sweep()

            if len(measured_power) != len(self.target_spectrum):
                print(f"警告: 测量光谱(长度{len(measured_power)})与目标光谱(长度{len(self.target_spectrum)})长度不匹配！")
                return 1e6

            cost = np.mean((measured_power - self.target_spectrum)**2)
            print(f"评估完成, 成本值为: {cost:.6f}")
            return cost
        except Exception as e:
            print(f"评估过程中发生错误: {e}")
            return 1e6
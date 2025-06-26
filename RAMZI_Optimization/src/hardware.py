# src/hardware.py
import numpy as np
import time
import pyvisa 

class VoltageSource:
    """控制电压源的类。"""
    def __init__(self, address: str, mock_mode: bool = False):
        self.mock_mode = mock_mode
        self.address = address
        self.instrument = None
        
        if not self.mock_mode:
            # rm = pyvisa.ResourceManager()
            # self.instrument = rm.open_resource(self.address)
            # self.instrument.timeout = 5000
            print(f"电压源已连接到: {self.address}")
        else:
            print("电压源以模拟模式初始化。")

    def set_all_voltages(self, voltages: list[float]):
        """设置所有通道的电压。"""
        if not self.mock_mode and self.instrument:
            for i, v in enumerate(voltages):
                # 伪代码，指令需要您根据手册来填充
                # self.instrument.write(f'SET_VOLT CHAN{i+1}, VOLT{v}')
                pass
        
        print(f"设置电压: {[f'{v:.3f}' for v in voltages]}")

class OSA:
    """控制光谱分析仪的类。"""
    def __init__(self, address: str, mock_mode: bool = False):
        self.mock_mode = mock_mode
        self.address = address
        self.instrument = None
        
        if not self.mock_mode:
            # rm = pyvisa.ResourceManager()
            # self.instrument = rm.open_resource(self.address)
            # self.instrument.timeout = 20000
            print(f"OSA已连接到: {self.address}")
        else:
            print("OSA以模拟模式初始化。")

    def _parse_data_string(self, data_str: str) -> np.ndarray:
        """将逗号分隔的字符串解析为Numpy数组。"""
        return np.array([float(p) for p in data_str.strip().split(',')])

    def sweep(self) -> tuple[np.ndarray, np.ndarray]:
        """执行一次波长扫描并返回 (波长数组, 功率数组)。"""
        from . import config

        if not self.mock_mode and self.instrument:
            self.instrument.write(f':sens:wav:star {config.OSA_WAVELENGTH_START}nm')
            # ... 其他真实硬件指令 ...
            self.instrument.write(':init')
            self.instrument.query('*OPC?')
            lamb_str = self.instrument.query(':trac:x? tra')
            power_str = self.instrument.query(':trac:y? tra')
            wavelengths = self._parse_data_string(lamb_str)
            powers = self._parse_data_string(power_str)
        else:
            # 模拟硬件行为
            print("模拟OSA扫描...")
            num_points = int((config.OSA_WAVELENGTH_STOP - config.OSA_WAVELENGTH_START) / config.OSA_RESOLUTION) + 1
            wavelengths = np.linspace(config.OSA_WAVELENGTH_START, config.OSA_WAVELENGTH_STOP, num_points)
            noise = np.random.normal(0, 0.5, num_points)
            peak_shift = np.random.uniform(-0.2, 0.2)
            powers = -25 * np.exp(-(wavelengths - (1550.5 + peak_shift))**2 / 0.1) - 10 + noise
            time.sleep(0.5)

        print("扫描完成。")
        return wavelengths, powers
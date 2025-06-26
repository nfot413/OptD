# tests/test_core_logic.py
import numpy as np
import pytest

# 导入src目录下的模块用于测试
from src.objective import ObjectiveFunction
from src.hardware import OSA, VoltageSource

def test_objective_function_with_mocks(mocker):
    """
    测试目标函数的核心逻辑，并模拟(mock)所有硬件。
    'mocker'是pytest-mock插件提供的工具。
    """
    # 1. 准备 (Arrange)
    # 创建硬件类的模拟实例
    mock_osa = mocker.Mock(spec=OSA)
    mock_vs = mocker.Mock(spec=VoltageSource)
    
    # “教”模拟的OSA在被调用.sweep()时返回什么
    fake_powers = np.full(10, -15.0)
    mock_osa.sweep.return_value = (None, fake_powers)
    
    # 准备目标光谱和输入电压
    target_spectrum = np.full(10, -17.0)
    test_voltages = [0.5] * 11
    
    # 用模拟硬件实例化目标函数
    objective = ObjectiveFunction(mock_osa, mock_vs, target_spectrum)
    
    # 2. 执行 (Act)
    cost = objective(test_voltages)
    
    # 3. 断言 (Assert)
    # 验证 set_all_voltages 方法是否被用正确的参数调用了一次
    mock_vs.set_all_voltages.assert_called_once_with(test_voltages)
    
    # 验证 osa.sweep 方法是否被调用了一次
    mock_osa.sweep.assert_called_once()
    
    # 验证成本计算是否正确: mean((-15.0 - (-17.0))^2) = mean(2^2) = 4.0
    assert cost == pytest.approx(4.0)
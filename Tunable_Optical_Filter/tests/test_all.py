import pytest
import numpy as np
from numpy.testing import assert_allclose

from src.components import TunableMZI, PhaseShifter, MicroRingResonator, DelayLine
from src.simulation import run_simulation, create_reference_box_filter

# ==============================================================================
# 1. 测试基础光学组件 (来自 components.py)
# ==============================================================================

def test_tunable_mzi_matrices():
    """测试MZI在特定相移下的传输矩阵是否符合物理预期。"""
    j = 1j
    coupler_50_50 = 0.5 * np.array([[-1+j, 1+j], [1+j, -1+j]])

    # Case 1: 当 theta = 0, MZI 处于“交叉”状态
    mzi_cross = TunableMZI(theta=0)
    # 此时内部相移矩阵是单位矩阵 I，所以 T = C @ I @ C = C @ C
    expected_cross_matrix = coupler_50_50 @ coupler_50_50
    # 经过计算，C @ C = [[-j, 0], [0, -j]]
    assert_allclose(mzi_cross.get_transfer_matrix(), expected_cross_matrix, atol=1e-9)

    # Case 2: 当 theta = pi, MZI 处于“直通”状态
    mzi_bar = TunableMZI(theta=np.pi)
    # 此时内部相移矩阵 P = [[-1, 0], [0, 1]]
    phase_matrix_pi = np.array([[-1, 0], [0, 1]])
    expected_bar_matrix = coupler_50_50 @ phase_matrix_pi @ coupler_50_50
    # 经过计算，结果为 [[j, 0], [0, -j]]，这是一个对角矩阵（直通状态）
    assert_allclose(mzi_bar.get_transfer_matrix(), expected_bar_matrix, atol=1e-9)

def test_phase_shifter_identity():
    """测试零相移的PhaseShifter是否为单位矩阵。"""
    ps_identity = PhaseShifter(phi_t=0, phi_b=0)
    assert_allclose(ps_identity.get_transfer_matrix(), np.identity(2))

def test_micro_ring_resonator_on_resonance():
    """测试MRR在共振点上的响应。"""
    t, k = 0.98, 0.1
    mrr = MicroRingResonator(t=t, k=k, phi_offset=np.pi)
    
    # 选择一个共振频率, 此时 2*w + phi_offset 是 2*pi 的整数倍
    # 例如, 2*w + pi = 2*pi  => 2w = pi => w = pi/2
    w_res = np.pi / 2
    
    # 此时 exp(-j*(2w+phi_offset)) = exp(-j*2pi) = 1
    expected_response = (np.sqrt(1 - k) - t**2) / (1 - t**2 * np.sqrt(1 - k))
    assert_allclose(mrr.get_response(w_res), expected_response)

def test_delay_line_at_zero_frequency():
    """测试延迟线在零频率(w=0)下的响应。"""
    dl = DelayLine(t=0.9, delay=1.0, phi_c=np.pi/2)
    # 当 w=0, 响应应为 t * exp(-j*phi_c)
    expected_response = 0.9 * np.exp(-1j * np.pi/2)
    assert_allclose(dl.get_response(w=0), expected_response)


# ==============================================================================
# 2. 测试仿真和辅助函数 (来自 simulation.py)
# ==============================================================================

def test_create_reference_box_filter():
    """测试理想方波滤波器的生成逻辑。"""
    freqs = np.linspace(100, 200, 101) # 从100到200共101个点
    pass_db, stop_db = 0, -40
    
    ref = create_reference_box_filter(
        frequency_array=freqs,
        center_freq=150,
        fsr=100, # 在此频率范围内无周期性
        bandwidth=20, # 通带应为 [140, 160]
        passband_level_db=pass_db,
        stopband_level_db=stop_db
    )
    
    assert ref.shape == freqs.shape
    # 中心点 (150) 应该在通带内
    assert ref[50] == pytest.approx(pass_db)
    # 边界点 (140, 160) 应该在通带内
    assert ref[40] == pytest.approx(pass_db)
    assert ref[60] == pytest.approx(pass_db)
    # 边界之外的点 (139, 161) 应该在阻带内
    assert ref[39] == pytest.approx(stop_db)
    assert ref[61] == pytest.approx(stop_db)

def test_run_simulation_simple_case():
    """使用一个极简的、可手动计算的场景测试核心仿真函数 run_simulation。"""
    # 场景: 上臂1个环, 下臂0个环。所有相移为0。
    n_ku, m_kl = 1, 0
    t = 0.98

    # 创建最简单的组件
    H1 = TunableMZI(theta=0).get_transfer_matrix()
    H3 = TunableMZI(theta=0).get_transfer_matrix()
    delay_line = DelayLine(t=t, delay=0, phi_c=0)
    upper_arm_mrrs = [MicroRingResonator(t=t, k=0, phi_offset=np.pi)]
    lower_arm_mrrs = []

    # 简单的频率范围
    w_range = np.array([0, np.pi])

    # 待“优化”的参数 (这里 k=0)
    params = [0.0] 

    # 运行仿真
    result_db = run_simulation(
        params=params, w_range=w_range, H1=H1, H3=H3,
        upper_arm_mrrs=upper_arm_mrrs, lower_arm_mrrs=lower_arm_mrrs,
        delay_line=delay_line, n_ku=n_ku, m_kl=m_kl
    )

    # 检查输出的基本属性
    assert isinstance(result_db, np.ndarray)
    assert result_db.shape == w_range.shape
    assert np.isreal(result_db).all(), "输出的dB值应为实数"

    # 根据正确的物理推导，H11的幅值应该是t
    expected_amplitude = t
    expected_db = 20 * np.log10(np.abs(expected_amplitude))

    # 用正确的期望值进行断言
    assert result_db[0] == pytest.approx(expected_db)
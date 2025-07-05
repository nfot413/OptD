import numpy as np

def run_simulation(params, w_range, H1, H3, upper_arm_mrrs, lower_arm_mrrs, delay_line, n_ku, m_kl):
    """
    执行光学滤波器的核心仿真计算。

    参数:
        params (list): 待优化的耦合系数k。
        w_range (np.ndarray): 角频率数组。
        H1, H3 (np.ndarray): 输入和输出部分的静态传输矩阵。
        upper_arm_mrrs, lower_arm_mrrs (list): 包含MRR实例的列表。
        delay_line (DelayLine): 延迟线实例。
        n_ku, m_kl (int): 上下臂的MRR数量。
    返回:
        np.ndarray: 仿真的功率传输谱 (dB)。
    """
    ku_params = params[:n_ku]
    kl_params = params[n_ku:]
    len_w = len(w_range)

    # 更新上臂MRR的耦合系数并计算响应
    for mrr, k in zip(upper_arm_mrrs, ku_params):
        mrr.update_coupling_coefficient(k)
    Au_mrr_responses = [mrr.get_response(w_range) for mrr in upper_arm_mrrs]
    Au = np.prod(Au_mrr_responses, axis=0)
    
    # 更新下臂MRR的耦合系数并计算响应
    if m_kl > 0:
        for mrr, k in zip(lower_arm_mrrs, kl_params):
            mrr.update_coupling_coefficient(k)
        Al_mrr_responses = [mrr.get_response(w_range) for mrr in lower_arm_mrrs]
        Al_mrr_product = np.prod(Al_mrr_responses, axis=0)
    else:
        Al_mrr_product = 1.0

    # 组合下臂的MRR和延迟线响应
    Al = Al_mrr_product * delay_line.get_response(w_range)

    # 构建中间矩阵 H2
    H2_stack = np.zeros((len_w, 2, 2), dtype=complex)
    H2_stack[:, 0, 0] = Au
    H2_stack[:, 1, 1] = Al

    # 计算最终传输函数
    H_final = H1 @ H2_stack @ H3
    H11 = H_final[:, 0, 0]
    
    return 20 * np.log10(np.abs(H11))

def objective_function(params, w_range, target_spectrum_db, H1, H3, upper_arm_mrrs, lower_arm_mrrs, delay_line, n_ku, m_kl):
    """
    优化器的目标函数（均方误差）。
    """
    simulated_spectrum = run_simulation(params, w_range, H1, H3, upper_arm_mrrs, lower_arm_mrrs, delay_line, n_ku, m_kl)
    mse_loss = np.mean((simulated_spectrum - target_spectrum_db)**2)
    return mse_loss

def create_reference_box_filter(frequency_array, center_freq, fsr, bandwidth, passband_level_db, stopband_level_db):
    """创建一个周期性的方形滤波器光谱作为目标。"""
    reference_signal = np.full_like(frequency_array, stopband_level_db)
    f_offset = frequency_array - center_freq
    f_normalized = np.mod(f_offset + fsr / 2, fsr) - fsr / 2
    passband_mask = np.abs(f_normalized) <= (bandwidth / 2)
    reference_signal[passband_mask] = passband_level_db
    return reference_signal
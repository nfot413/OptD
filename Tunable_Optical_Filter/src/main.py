import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
from functools import partial

# 从同级模块中导入：配置、组件类、仿真函数
from . import config
from .components import TunableMZI, PhaseShifter, MicroRingResonator, DelayLine
from .simulation import run_simulation, create_reference_box_filter, objective_function

def setup_components():
    """根据配置创建并返回所有光学组件的实例。"""
    # 预计算静态矩阵 H1 和 H3
    mzi_in = TunableMZI(config.THETA_I)
    mzi_out = TunableMZI(config.THETA_O)
    phase_shifter = PhaseShifter(config.PHI_T, config.PHI_B)
    H1 = mzi_in.get_transfer_matrix() @ phase_shifter.get_transfer_matrix()
    H3 = mzi_out.get_transfer_matrix()

    # 创建动态组件（即需要参与循环的组件）
    delay_line = DelayLine(config.T_SELF_COUPLING, config.DELAY_VAL, config.DELAY_PHASE)
    upper_arm_mrrs = [MicroRingResonator(t=config.T_SELF_COUPLING, k=0, phi_offset=np.pi) for _ in range(config.N_KU)]
    lower_arm_mrrs = [MicroRingResonator(t=config.T_SELF_COUPLING, k=0, phi_offset=np.pi) for _ in range(config.M_KL)]

    # 将所有组件打包到一个字典中，方便传递
    components = {
        'H1': H1, 'H3': H3, 'delay_line': delay_line,
        'upper_arm_mrrs': upper_arm_mrrs, 'lower_arm_mrrs': lower_arm_mrrs,
    }
    return components

def prepare_environment():
    """准备仿真的频率数组和目标频谱。"""
    # 频率数组
    w1 = config.W_MIN_PI_MULTIPLE * np.pi
    w2 = config.W_MAX_PI_MULTIPLE * np.pi
    w_range = np.arange(w1, w2, config.DW_STEP)
    len_w = len(w_range)
    frequency_f = np.linspace(
        config.F_CENTER - config.S_FSR_SPAN * config.FSR,
        config.F_CENTER + config.S_FSR_SPAN * config.FSR,
        len_w
    )

    # 目标频谱
    target_spectrum = create_reference_box_filter(
        frequency_array=frequency_f,
        center_freq=config.F_CENTER,
        fsr=config.FSR,
        bandwidth=config.BANDWIDTH,
        passband_level_db=config.PASSBAND_DB,
        stopband_level_db=config.STOPBAND_DB
    )

    # 将环境参数打包到字典中
    environment = {
        'w_range': w_range,
        'frequency_f': frequency_f,
        'target_spectrum': target_spectrum
    }
    return environment

def run_and_display_results(components, environment):
    """执行优化，并打印、绘制最终结果。"""
    w_range = environment['w_range']
    target_spectrum = environment['target_spectrum']
    
    # 使用 functools.partial 创建一个更简单的目标函数。
    # 新函数 objective_for_optimizer 只有一个参数 `params`，这正是优化器所需要的。
    # objective_function 的其他参数都被“冻结”或“预设”了。
    objective_for_optimizer = partial(
        objective_function,
        w_range=w_range,
        target_spectrum_db=target_spectrum,
        H1=components['H1'], H3=components['H3'],
        upper_arm_mrrs=components['upper_arm_mrrs'],
        lower_arm_mrrs=components['lower_arm_mrrs'],
        delay_line=components['delay_line'],
        n_ku=config.N_KU, m_kl=config.M_KL
    )

    # --- 运行优化 ---
    total_params = config.N_KU + config.M_KL
    bounds = [(0, 1)] * total_params
    print("="*50)
    print(f"开始优化 {total_params} 个参数 ({config.N_KU}个上臂环, {config.M_KL}个下臂环)...")
    print("="*50)
    
    start_time = time.time()
    result = differential_evolution(
        objective_for_optimizer,
        bounds,
        maxiter=config.MAX_ITERATIONS,
        popsize=config.POPULATION_SIZE,
        disp=True,
        workers=-1
    )
    end_time = time.time()
    print(f"\n优化完成！总耗时: {end_time - start_time:.2f} 秒")
    
    # --- 显示结果 ---
    best_params = result.x
    print("\n" + "="*50)
    print("优化结果详情:")
    print(f"  - 最低损失值 (MSE): {result.fun:.6f}")
    print(f"  - 找到的最佳参数 ({total_params}个):")
    for i in range(config.N_KU):
        print(f"    ku{i+1} = {best_params[i]:.4f}")
    for i in range(config.M_KL):
        print(f"    kl{i+1} = {best_params[i + config.N_KU]:.4f}")
    print("="*50)
    
    # --- 最终验证与绘图 ---
    print("\n正在使用找到的最佳参数进行最终效果验证...")
    final_spectrum_db = run_simulation(
        params=best_params, w_range=w_range, **components, n_ku=config.N_KU, m_kl=config.M_KL
    )

    plt.figure(figsize=(14, 7))
    plt.plot(environment['frequency_f'], target_spectrum, 'r--', lw=2.5, label='理想方波目标')
    plt.plot(environment['frequency_f'], final_spectrum_db, 'b-', lw=1.5, label='优化后参数得到的最终响应')
    plt.title(config.PLOT_TITLE, fontsize=16)
    plt.xlabel(config.PLOT_XLABEL, fontsize=12)
    plt.ylabel(config.PLOT_YLABEL, fontsize=12)
    plt.grid(True)
    plt.ylim(config.PLOT_YLIM_MIN, config.PLOT_YLIM_MAX)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()

def main():
    """项目主流程控制器。"""
    # 配置matplotlib以支持中文显示
    plt.rcParams['font.sans-serif'] = ['Songti SC']
    plt.rcParams['axes.unicode_minus'] = False

    # 1. 创建光学组件
    components = setup_components()
    
    # 2. 准备仿真环境
    environment = prepare_environment()
    
    # 3. 执行优化并显示结果
    run_and_display_results(components, environment)

if __name__ == '__main__':
    # 此脚本应作为模块(`python -m src.main`)运行，以确保相对导入正常工作
    main()
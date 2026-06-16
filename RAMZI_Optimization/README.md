# RAMZI 仿真与优化

本项目用于 RAMZI 架构建模、光谱计算与耦合系数优化。

## 功能

- 构建包含方向耦合器、全通环形谐振器和下臂延迟线的 RAMZI 传输函数模型。
- 从 `data/fsr1.xlsx` 读取目标光谱。
- 计算 RAMZI 输出端口的仿真光谱，单位为 dB。
- 使用 SciPy 差分进化算法优化环形谐振器耦合系数。
- 保存优化结果数据和光谱对比图。

## 目录结构

```text
RAMZI_Optimization/
  README.md
  requirements.txt
  pytest.ini
  data/
    fsr1.xlsx
  src/
    __init__.py      # 包导出入口
    components.py    # 光学组件传输函数
    config.py        # 模型、数据、优化器和输出配置
    experiment.py    # 仿真优化实验流程
    io.py            # 目标光谱读取和结果保存
    main.py          # 命令行入口
    objective.py     # 损失函数和优化目标
    optimizer.py     # 优化器辅助函数
    plotting.py      # 绘图输出
    simulation.py    # RAMZI 传输函数仿真
  tests/
    test_components.py    # 光学组件测试
    test_experiment.py    # 实验流程测试
    test_simulation.py    # 仿真与目标函数测试
```

## 主要模块

- `components.py` 定义方向耦合器、全通环形谐振器和延迟线的基础传输函数。
- `simulation.py` 根据 RAMZI 结构参数计算复数传输函数和 dB 光谱。
- `objective.py` 定义均方误差损失和可调用优化目标。
- `optimizer.py` 封装差分进化优化流程。
- `experiment.py` 串联频率网格、目标光谱、优化器、结果保存和绘图。
- `config.py` 提供数据、模型、频率网格、优化器和输出配置。

## 依赖

依赖列表位于 `requirements.txt`：

```text
numpy
pandas
openpyxl
scipy
matplotlib
pytest
```

安装依赖：

```bash
pip install -r requirements.txt
```

## 运行

在 `RAMZI_Optimization` 目录执行：

```bash
python -m src.main
```

输出文件：

```text
results/optimized_spectrum.csv
results/optimized_spectrum.png
```

## 测试

在 `RAMZI_Optimization` 目录执行：

```bash
pytest -q
```

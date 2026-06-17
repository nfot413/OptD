# RAMZI 仿真与优化

本项目用于固定 RAMZI 结构的传输函数仿真、方波目标构造与相位参数优化。

## 功能

- 构建参考 MATLAB 原型的固定 RAMZI 结构。
- 使用输入 3dB MZI、四个微环单元、上下臂相位和输出 3dB MZI 计算输出光谱。
- 根据微环群折射率和环长计算每个微环的 FSR。
- 生成中心频率附近的方波目标光谱。
- 优化 `fai1`、`fai2`、`fai3`、`fai4`、`fait`、`theta13`、`theta24` 七个相位参数。
- 保存优化结果数据和光谱对比图。
- 提供本地图形界面，用于修改参数、运行仿真和执行优化。

## 固定结构

RAMZI 结构包含：

- 输入 MZI：`Ki = 0.5`，`thetai = pi/2`。
- 上臂：`Ring 1`、`Ring 2` 和整体相位 `fait`。
- 下臂：`Ring 3`、`Ring 4` 和整体相位 `faib = -fait`。
- 输出 MZI：`Ko = 0.5`，`thetao = pi/2`。

微环物理参数：

```text
Ring 1: L1 = 350 um,  ng = 4.3
Ring 2: L2 = 3000 um, ng = 4.3
Ring 3: L3 = 350 um,  ng = 4.3
Ring 4: L4 = 3000 um, ng = 4.3
```

FSR 使用：

```text
FSR = c / (ng * L)
```

其中 `c = 3e8 m/s`，`L` 为微环一圈的等效长度。

## 优化参数

优化变量为：

```text
fai1
fai2
fai3
fai4
fait
theta13  # theta1 = theta3
theta24  # theta2 = theta4
```

约束关系：

```text
faib = -fait
theta1 = theta3 = theta13
theta2 = theta4 = theta24
```

参数范围：

```text
fai1, fai2, fai3, fai4, fait, theta13, theta24 ∈ [-pi, pi]
```

## 目录结构

```text
RAMZI_Optimization/
  README.md
  requirements.txt
  pytest.ini
  gui/
    app.js
    index.html
    styles.css
  src/
    architecture.py  # 固定 RAMZI 结构、物理参数、方波目标和优化入口
    __init__.py      # 包导出入口
    components.py    # 基础光学组件函数
    config.py        # 默认频率、优化器和输出配置
    experiment.py    # 命令行仿真优化流程
    gui_server.py    # 本地图形界面服务
    io.py            # 结果保存
    main.py          # 命令行入口
    objective.py     # 损失函数
    optimizer.py     # 优化器辅助函数
    plotting.py      # 绘图输出
    simulation.py    # 基础 RAMZI 传输函数工具
  tests/
    test_architecture.py  # 固定结构测试
    test_components.py    # 光学组件测试
    test_experiment.py    # 实验流程测试
    test_gui_server.py    # 图形界面接口测试
    test_simulation.py    # 基础仿真测试
```

## 依赖

依赖列表位于 `requirements.txt`：

```text
numpy
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

## 图形界面

在 `RAMZI_Optimization` 目录执行：

```bash
python -m src.gui_server
```

运行后会自动打开：

```text
http://127.0.0.1:8765
```

只启动服务、不自动打开网页：

```bash
python -m src.gui_server --no-open
```

界面支持：

- 显示固定 RAMZI 结构示意图。
- 显示每个微环的长度、群折射率和 FSR。
- 修改七个待优化相位参数。
- 设置中心频率、扫描范围、方波通带宽度和采样点数。
- 显示方波目标光谱与 RAMZI 仿真光谱。
- 执行参数优化并更新参数。

## 测试

在 `RAMZI_Optimization` 目录执行：

```bash
pytest -q
```

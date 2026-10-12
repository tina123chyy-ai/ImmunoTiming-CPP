# 零基础操作指南：如何在 VS Code 里运行

**这不是一个网页项目**；它是一个可运行的 C++17 命令行程序。项目里的Python脚本只负责把C++生成的CSV画成图。

## 你要准备什么

- 安装 VS Code，以及扩展 **C/C++**（Microsoft）与 **CMake Tools**。
- macOS：先安装 Xcode Command Line Tools（终端运行 `xcode-select --install`），并安装 CMake。
- Windows：安装 Visual Studio Build Tools（勾选 C++工具）和 CMake；推荐使用 Visual Studio 的 Developer PowerShell。
- 可选：如果你想生成PNG图，再安装 Python 和 `matplotlib`。

## 第一次运行（macOS/Linux）

1. 解压整个 `ImmunoTiming-CPP.zip`。
2. 在VS Code选择 **File → Open Folder**，选择 `ImmunoTiming-CPP` 文件夹。
3. 选择菜单 **Terminal → New Terminal**（这个终端属于VS Code，不是外部软件）。
4. 复制粘贴以下命令，每一行按回车：

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/immunotiming_cli --mode baseline --out results/my_first_run
```

5. 左边文件栏打开 `results/my_first_run/summary.csv`，查看四种方案的示意性计算结果。
6. 如已安装Matplotlib，可以生成图：

```bash
python3 -m pip install matplotlib
python3 visualization/plot_results.py results/my_first_run --output results/my_first_run.png
```

## 我应该先看哪些C++代码？

1. `include/TumorModel.h`：学习 `struct`、参数和状态变量。
2. `src/TumorModel.cpp`：学习如何把3条微分方程写成C++。
3. `src/ODESolver.cpp`：学习 `RK4` 的数值计算原理。
4. `src/TreatmentScheduler.cpp`：学习优先队列与事件调度。
5. `src/ScheduleRanker.cpp`：学习堆实现Top-K筛选。

## 亲自完成的任务（很重要）

如果你希望把项目作为自己的编程能力证明，**不能只是将现成代码上传GitHub**。建议先独立复写`TumorModel.cpp`中的三个导数公式，再尝试独立实现一个Euler求解器和一个简单的优先队列事件调度器。之后增加你自己写的测试和实验说明，在Git commit记录里保留这些工作。

## 常见问题

- `cmake: command not found`：还没安装CMake或者没有加入PATH。
- `No CMAKE_CXX_COMPILER could be found`：缺少C++编译器。
- 绘图报`No module named matplotlib`：安装对应Python环境中的matplotlib。
- 终端里显示了数值但没有图：C++默认输出CSV；需要可选Python脚本绘图。

**特别提醒：** 项目所有模型参数均为虚构且归一化的演示值。结果只可用于学习数学建模和程序设计，不能用于真实患者治疗计划。

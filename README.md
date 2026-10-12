# ImmunoTiming-CPP 🧪
**A C++17 simulation and treatment-event scheduling framework for a conceptual CAR-T + personalized neoantigen vaccine model**

> **Research/education only.** All parameter values, rates, and results are normalized **hypothetical examples**. No patient data, clinical parameter fitting, safety model, or treatment recommendation is included. The "best schedule" is the best *under the toy equations and a deliberately arbitrary objective*, not a treatment schedule for people.

## 项目简介 / Overview

这个项目延伸一篇关于**个性化新抗原疫苗 + CAR-T 联合治疗实体瘤**的综述与概念性数学模型。原研究使用三个状态变量描述肿瘤负荷、CAR-T细胞和疫苗诱导的效应T细胞，辅以假设数据与示意曲线。此项目新增的是**可运行的 C++ 数值求解器、治疗事件调度算法、Top-K 堆排序和参数敏感性示例**。

The original paper's conceptual equations are implemented here, but **the hypothetical parameters in this repository are newly selected demonstration values**, not values measured in patients. The program does not establish biological validity.

### Research question / 研究问题
How do different hypothetical administration times change the behavior of a simplified, coupled differential-equation model?

研究对象是数学模型中的时间安排对**模型结果**的影响，而不是临床最佳治疗顺序。

## 功能 / Features

- **C++17** 3-state ODE model, numerical Euler and 4th-order Runge–Kutta (RK4) solvers
- `std::priority_queue` treatment-event scheduler (start/end of finite-rate pulses)
- Four baseline scenarios: no treatment, vaccine alone, CAR-T alone, combined
- Candidate schedule ranking using a bounded **Top-K max-heap**
- Hypothetical `beta` (antigen-release feedback) and `gamma` (vaccine-to-CAR-T synergy) sensitivity grid
- CSV exports + optional Python/Matplotlib plotter
- CMake/CTest automatic checks

## 快速运行 / Quick start

You need a C++17 compiler (GCC/Clang/MSVC) and CMake 3.16+.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/immunotiming_cli --mode baseline --out results/baseline
./build/immunotiming_cli --mode rank --k 5 --out results/rank
./build/immunotiming_cli --mode sensitivity --out results/sensitivity
```

Windows PowerShell (if using Visual Studio):

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
.\build\Release\immunotiming_cli.exe --mode baseline --out results\baseline
```

On macOS, after installing Xcode Command Line Tools and CMake, the Unix commands above work.

Plot figures (optional, `pip install matplotlib` required):

```bash
python3 visualization/plot_results.py results/baseline --output results/baseline.png
```

### 命令行示例 / CLI

```bash
./build/immunotiming_cli --help
./build/immunotiming_cli --mode single --vaccine-day 2 --cart-day 12 --out results/custom
./build/immunotiming_cli --mode baseline --method euler --step 0.01 --days 40 --out results/euler
```

Available modes: `baseline`, `single`, `rank`, `sensitivity`.

- `--vaccine-day` and `--cart-day`: *model-time* start for each hypothetical input
- `--days`: model-time horizon (default 40)
- `--step`: integration step (default 0.05)
- `--method`: `rk4` (default) or `euler`
- `--k`: number of retained schedules in ranking (default 5)

## 核心方程 / Original conceptual model

Denote tumor burden by **T**, CAR-T population by **C**, and endogenous effector population by **E**.

```text
dT/dt = r*T*(1 - T/K) - kC*C*T - kE*E*T

dC/dt = alphaC*(1+gamma*E)*(T/(hC+T))*C - deltaC*C + uC(t)

dE/dt = alphaV*V(t) + beta*C*T - deltaE*E
```

- `gamma`: hypothesized **E→C** positive interaction in CAR-T expansion
- `beta`: hypothesized **C×T→E** feedback from CAR-T-mediated tumor killing/antigen release
- `V(t)`: vaccine input rate; `uC(t)`: CAR-T input rate
- Inputs are **finite-duration rectangular pulses** for computational convenience: 2 model-time units for vaccine and 1 for CAR-T. They are **not** medical dosing protocols.

The paper's Figure 2 is illustrative, manually drawn from hypothetical concepts. **This repository's curves are newly computed ODE trajectories**, and they need not reproduce that figure. See [`docs/MODEL_AND_LIMITATIONS.md`](docs/MODEL_AND_LIMITATIONS.md).

## 数据结构设计 / Data structures

| Structure | Real use | Complexity |
|---|---|---|
| `std::vector<Point>` | Sampled trajectory | amortized O(1) append |
| `std::priority_queue<Event>` | Chronological treatment start/end events | O(log m) per event push/pop |
| `std::priority_queue<Candidate>` | Keep best K candidate schedules | O(log K) per candidate after evaluation |
| `std::sort` | Final ranked results | O(K log K) |

For N schedule candidates, ranking overhead is `O(N log K)` **plus O(N × S)** ODE simulation work (S is the number of solver steps per candidate). In these experiments, numerical integration dominates runtime.

## 仓库结构 / Repository

```text
include/         C++ public API and model/data types
src/             Model, RK4/Euler solver, treatment scheduling, ranking, CSV export, CLI
tests/           Unit/behavior tests; CTest
visualization/   Optional Python CSV plotting (no biological computation)
docs/            Model assumptions, limitations, technical design, study guide
examples/        Example command recipes
results/         Example generated CSV and figure(s)
```

## 如何读懂和改进 / Student learning path

1. Read `include/TumorModel.h` then `src/TumorModel.cpp` to understand the equations.
2. Read `src/ODESolver.cpp`, derive RK4's 4 intermediate slopes, compare with Euler.
3. Read `src/TreatmentScheduler.cpp`, explain why the program *splits integration steps at event boundaries*.
4. Read `src/ScheduleRanker.cpp`, derive the Top-K heap time complexity.
5. Try `--mode sensitivity`, and explain why changing `beta`/`gamma` only changes **hypothetical model behavior**.
6. Reproduce a plot and modify the objective function; document how the ranking changes.

## Attribution & academic integrity / 署名与学术诚信

The accompanying conceptual equations are adapted from the student's source manuscript, *Personalized Neoantigen Vaccine Combined with CAR-T Therapy for Solid Tumor Immunotherapy*. This repository is a **new educational programming implementation** derived from those equations, not code or simulation produced during the original manuscript's writing. When presenting the work, distinguish the original modeling contribution from any assistance received for the C++ code, and claim only parts you independently implemented, debugged, or understood.

## License
MIT. See `LICENSE`.

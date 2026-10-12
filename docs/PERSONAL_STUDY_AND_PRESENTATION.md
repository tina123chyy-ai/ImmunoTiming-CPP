# 学习任务与申请展示建议

## 你已经完成的旧研究 vs 本项目新完成的工作

原综述论文：独立选题、以生物学机制构建概念性微分方程模型、使用假设参数和手绘示意曲线。它未进行患者数据拟合，也不是临床疗效验证。

本GitHub项目：将同一结构的数学方程编码进C++17，计算RK4轨迹、管理治疗事件、比较示意性时间安排并进行Top-K与敏感性示例分析。属于**后续编程实现**。

## 不能直接声称的内容

- “使用AI预测患者最佳治疗时间”——这里没有机器学习和真实临床预测。
- “用模拟验证CAR-T联合疫苗有效”——演示结果取决于人为设定的正向协同项。
- “原论文已经使用C++求解微分方程”——代码是随后新制作的。
- “所有代码由我独立实现”——除非你确实自行完成或在辅导后实质性重写并能解释。

## 建议你亲自完成的四项任务

1. 在纸上推导RK4的4个斜率，并在`ODESolver.cpp`逐行备注解释。
2. 不看源码，自己重写`TreatmentScheduler`的最小版本并用两种给药事件测试。
3. 向老师解释为什么Top-K要维护一个“堆顶是最差者”的最大堆，而不是每次全部排序。
4. 亲自尝试修改 `beta` 为0，写下对结果影响的分析与原因。

## 英文简历/活动栏内容（需真实完成后再使用）

Developed a C++17 implementation of a conceptual three-state immunotherapy ODE model, including RK4 integration, heap-based treatment-event scheduling, and Top-K ranking of hypothetical treatment scenarios; documented sensitivity to assumed feedback terms and limits of clinical interpretation.

## GitHub README优秀标准

- 目的明确；原论文与后续代码贡献区分清楚
- 一条命令可完成构建，自动测试通过
- 实验流程可复现，图中注明假设参数
- 含数据结构选择理由和时间复杂度分析
- 公开模型局限，而不是只呈现看起来漂亮的联合治疗曲线

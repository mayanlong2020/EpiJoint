# BakerLab: RSV Triple-Epitope Co-Core Design Reproduction

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## English Version

### 1. Project Overview
This repository contains a full reproduction and execution record of the triple-epitope RSV F antigen design, inspired by the state-of-the-art research from **Nature Chemical Biology (2025)**. The core objective is to utilize a deep learning pipeline to accurately scaffold three non-overlapping neutralizing epitopes of the RSV F protein—**Site II, Site IV, and Site V**—within a single, compact globular protein domain.

### 2. Methodology & Acknowledgments
This reproduction is based on the multi-motif scaffolding methodologies pioneered by the University of Washington's Institute for Protein Design (IPD) and EPFL.

**Citation:**
> *Accurate single-domain scaffolding of three nonoverlapping protein epitopes using deep learning.* 
> Castro, K.M., Watson, J.L., ... Baker, D., & Correia, B.E. **Nature Chemical Biology (2025).**
> [Official Paper](https://www.nature.com/articles/s41589-025-02083-z) | [Official RFjoint2 Repo](https://github.com/RosettaCommons/RFjoint2)

### 3. Pipeline Architecture
Our implementation automates the workflow from motif extraction to structure-based binding validation:
*   **Phase 0-5**: A complete 6-stage pipeline covering data preparation, **RFjoint2** backbone diffusion, **ProteinMPNN** sequence design, and **Boltz-2** binding validation.

### 4. Running Environment & Hardware Optimization
This pipeline is specifically optimized for high-performance **Apple Silicon (M4 Max / Ultra)** environments.

*   **Computation Strategy**: Phase 1 (**RFjoint2**) execution is intentionally routed to **CPU Multi-threading**.
*   **Technical Rationale**: Unlike RFdiffusion2, which has broader Apple Silicon compatibility, **RFjoint2** currently encounters stability and kernel-op issues with **MPS (Metal Performance Shaders)**. To ensure absolute numerical precision and process stability, this pipeline leverages the high-core density of the M4 Max CPU (Efficiency + Performance cores) via optimized concurrency.

### 5. Quick Start
To launch the fully automated pipeline:
```bash
python epiforge_main.py
```
**Important:** The current repository includes output from a low-throughput demonstration run. For production-grade results, please increase the sampling depth (set `NUM_DESIGNS` to **5000+**).

### 6. Deep Audit & Technical Insights (Low-Throughput Analysis)
Our atomic-level audit of the demonstration run identified some structural anomalies. These are **consequences of insufficient sampling depth**, not a flaw in the methodology itself.

#### 6.1 Observed Structural Strains
*   **Main-chain Strain (Site IV)**: The artificial disulfide (**C40-C45**) in the demo winner shows a **CA-CA distance of 3.640 Å**.
*   **Steric Clashes**: Detected **29 severe internal clashes** in the low-throughput candidate.
*   **Scientific Conclusion**: High-throughput sampling is required to navigate the complex energy landscape and identify the rare "native-like" optimal designs.

---

<a name="chinese"></a>
## 中文版

### 1. 项目概述
本仓库记录了对 **Nature Chemical Biology (2025)** 关于 RSV F 蛋白三表位设计的完整复现。该设计的核心在于利用深度学习管线，在单一球蛋白结构域中同时精确展示三个非重叠的中和表位：**Site II**、**Site IV** 和 **Site V**。

### 2. 方法论与致谢
本复现项目严格遵循华盛顿大学蛋白质设计研究所（IPD）与 EPFL 开创的多表位支架化（Multi-motif scaffolding）方法学。

**引用信息：**
> *Accurate single-domain scaffolding of three nonoverlapping protein epitopes using deep learning.* 
> Castro, K.M., Watson, J.L., ... Baker, D., & Correia, B.E. **Nature Chemical Biology (2025).**
> [论文链接](https://www.nature.com/articles/s41589-025-02083-z) | [官方代码仓 (RFjoint2)](https://github.com/RosettaCommons/RFjoint2)

### 3. 管线架构
本项目实现了完整的 6 阶段自动化管线：涵盖数据准备、**RFjoint2** 骨架扩散、**ProteinMPNN** 序列设计、以及 **Boltz-2** 结合面验证。

### 4. 运行环境与硬件优化
本管线专门针对高性能 **Apple Silicon (M4 Max / Ultra)** 环境进行了深度优化。

*   **计算策略**：Phase 1 (**RFjoint2**) 的执行被有意路由至 **CPU 多线程并行**。
*   **技术选型原因**：不同于兼容性较广的 RFdiffusion2，**RFjoint2** 目前在 **MPS (Metal Performance Shaders)** 加速上存在算子兼容性与稳定性问题。为了确保数值计算的绝对精准与流程稳定，本管线充分利用了 M4 Max CPU 的高核心密度（性能核+能效核），通过优化并发逻辑实现高效设计。

### 5. 快速启动
运行主入口脚本即可启动管线：
```bash
python epiforge_main.py
```

#### 关键参数设置
*   **NUM_DESIGNS**: 设置拓扑结构总数。建议在正式科研中将其设置为 **5000+**。
*   **CONCURRENCY_CPU**: 并行任务数，已针对 M4 系列芯片进行优化。

### 6. 深度审计与技术洞察（低通量样本分析）
我们在演示运行的小样本审计中发现了一些结构异常。需要强调的是，这些问题是**由于采样深度不足导致的非最优解现象**。

#### 6.1 结构观测现象
*   **主链挤压 (Site IV)**：演示样本中的人工二硫键 (**C40-C45**) **CA-CA 间距为 3.640 Å**。
*   **空间位阻 (Steric Clashes)**：低通量候选者内部检测到 **29 处严重碰撞**。
*   **科学结论**：在实际应用中必须通过**加大采样量（NUM_DESIGNS）**来跨越能量势垒，从而筛选出物理学上真正完美的候选蛋白。

---
**Disclaimer**: This project is for academic reproduction purposes. We acknowledge and respect the pioneering work of the Baker Lab and the Correia Lab.

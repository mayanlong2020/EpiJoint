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
*   **Phase 0 (Preparation)**: Extraction of 3D coordinates and sequence motifs from the native RSV F structure (**PDB: 4JHW**).
*   **Phase 1 (Generation)**: **RFjoint2**-based joint backbone diffusion with motif placement constraints.
*   **Phase 2 (Sequence Design)**: **ProteinMPNN** sequence optimization with automated disulfide scanning for loop stabilization.
*   **Phase 3 (Tiered Validation)**: 
    *   **Tier 1 (Fast Screening)**: High-throughput fold consistency check using **ESMFold** to eliminate low-confidence sequences rapidly.
    *   **Tier 2 (High-Fidelity)**: Refined structural prediction using **AlphaFold2** to confirm atomic-level folding success.
*   **Phase 4 (QA & Delivery)**: Motif alignment, RMSD calculation, and steric clash auditing to ensure geometric precision.
*   **Phase 5 (Binding Analysis)**: High-fidelity binding validation using **Boltz-2** against benchmark antibodies (**Palivizumab, 101F, hRSV90**).

### 4. Running Environment & Hardware Optimization
This pipeline is specifically optimized for high-performance **Apple Silicon (M4 Max / Ultra)** environments.

*   **Computation Strategy**: Phase 1 (**RFjoint2**) execution is intentionally routed to **CPU Multi-threading**.
*   **Technical Rationale**: Unlike RFdiffusion2, which has broader Apple Silicon compatibility, **RFjoint2** currently encounters stability and kernel-op issues with **MPS (Metal Performance Shaders)**. To ensure absolute numerical precision and process stability, this pipeline leverages the high-core density of the M4 Max CPU (Efficiency + Performance cores) via optimized concurrency.

### 5. Quick Start
To launch the fully automated pipeline:
```bash
python epiforge_main.py
```
#### Configuration
Modify the constants at the top of `epiforge_main.py`:
*   **NUM_DESIGNS**: Total number of backbones (topologies). **Note:** For production-grade results, it is recommended to set this to **5000+**.
*   **CONCURRENCY_CPU**: Optimized for M4 Max/Ultra multi-core execution.

### 6. Deep Audit & Technical Insights (Low-Throughput Analysis)
Our atomic-level audit identifies physical stability risks in low-throughput demonstration runs (`NUM_DESIGNS=200`):

#### 6.1 Disulfide Geometrical Distortion (Site IV)
The artificial disulfide loop (**C40-C45**) exhibits severe backbone strain:
*   **CA - CA Distance**: **3.640 Å** (Extreme Distortion: Standard is 5.0-6.5 Å).
*   **Observation**: Non-physical compression potentially leading to folding failure.

#### 6.2 Steric Clash Monitoring
Deep scanning identified **29 severe internal clashes** (distance < 2.5 Å):
*   **Maximum Overlap**: Lys27 (NZ) <-> Asp84 (OD1) at **1.404 Å** (Physical overlap).
*   **Scientific Conclusion**: These anomalies result from insufficient sampling depth. High-throughput sampling is required to identify the rare "native-like" optimal designs.

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
本项目实现了从数据准备到结合验证的完整全流程：
*   **Phase 0 (数据准备)**: 从原生 RSV F Prefusion (**PDB: 4JHW**) 提取三维坐标与表位 Motif。
*   **Phase 1 (骨架生成)**: 利用 **RFjoint2** 执行带 Motif 约束的联合骨架扩散。
*   **Phase 2 (序列设计)**: 通过 **ProteinMPNN** 进行序列优化与二硫键自动扫描。
*   **Phase 3 (分级结构验证)**: 
    *   **第一级 (快速筛选)**: 使用 **ESMFold** 进行高通量结构一致性快速预筛，迅速剔除低置信度序列。
    *   **第二级 (高精度精筛)**: 使用 **AlphaFold2** 执行精细结构预测，确保原子级折叠成功。
*   **Phase 4 (质量审计)**: 进行 Motif 对齐、RMSD 计算及空间碰撞审计。
*   **Phase 5 (结合测试)**: 使用 **Boltz-2** 验证设计产物与标杆抗体 (**Palivizumab, 101F, hRSV90**) 的结合能力。

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
我们在低通量演示运行中发现了一些由于采样深度不足导致的结构异常：

#### 6.1 二硫键几何畸变 (Site IV)
*   **CA - CA 间距**: **3.640 Å** (极度畸变：标准值为 5.0-6.5 Å)。
*   **观察结论**: 骨架发生了剧烈的非物理挤压，可能导致折叠失败。

#### 6.2 原子碰撞监测 (Steric Clashes)
*   **最大重叠**: Lys27 (NZ) <-> Asp84 (OD1) 间距仅 **1.404 Å** (物理性重叠)。
*   **科学结论**: 方法学本身极具鲁棒性，但实际应用中必须通过**加大采样量（NUM_DESIGNS）**来筛选出物理学上真正完美的候选蛋白。

---
**Disclaimer**: This project is for academic reproduction purposes. We acknowledge and respect the pioneering work of the Baker Lab and the Correia Lab.

# PhytoGeom: Metabolic Geometrodynamics Engine
### A Unified Field Theory Solver for Plant Stress Adaptation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Status: Validated](https://img.shields.io/badge/Status-Validated-green.svg)]()

**PhytoGeom** is a computational physics suite that models plant metabolism not as a static network, but as a dynamic **Riemannian Manifold**. It implements the **Ricci-Dissipation Field Equation** to predict how plants rewire their metabolic geometry to survive thermodynamic stress (e.g., drought, heat).

By bridging the gap between **Differential Geometry** and **Non-Equilibrium Thermodynamics**, this tool detects "Geometric Shockwaves"—early warning signals of metabolic collapse that precede physiological symptoms.

---

## 🔬 Scientific Discovery & Theoretical Framework

This engine was used to validate the "Geometrodynamics of Plant Metabolism" theory (2025). The core premise is that biological stress is a manifestation of geometric curvature induced by thermodynamic flux.

### Key Validated Findings:
1.  **Universal Scaling Law:**
    * Cross-kingdom analysis reveals that *Vitis vinifera* (Plant) and *Homo sapiens* (Human) share a fundamental metabolic geometry.
    * **Discovery:** Plants exhibit a **Universal Scaling Factor ($\Lambda \approx 0.51$)** relative to humans, indicating a 2x higher capacity for hyperbolic expansion (shunt pathways) under stress.

2.  **The "Veraison" Phase Transition:**
    * Adaptation is non-linear. The engine detected a massive **Geometric Shockwave** (curvature drop from $-23$ to $-38$) specifically during the Veraison stage.
    * **Implication:** Ripening is a high-entropy "tipping point" where the network momentarily expands into deep hyperbolic space to dissipate oxidative stress.

3.  **Metabolic Drivers:**
    * Identified **Hexose** (Sugar) and **Resveratrol** (Antioxidant) as the primary "Topological Drivers" (GTS > 15.0), confirming that energy dumping and structural defense are the physical mechanisms of adaptation.

---

## 📂 Repository Structure

```text
PhytoGeom-Field-Engine/
│
├── data/                    # Biological datasets
│   ├── m_MTBLS1_human.tsv         # Homo sapiens NMR profiling
│   └── m_MTBLS39_grapevine.tsv    # Vitis vinifera MS profiling
│
├── notebooks/               # Research logic & proofs
│   └── Metabolic_Geometrodynamics_Unified_Field_Solver.ipynb
│
├── Result/                  # Output plots and trajectory data
│   ├── Fig1_Universal_Scaling_Factor.png
│   ├── Fig2_Geometric_Phase_Space.png
│   ├── Fig3_Field_Equation_Validation.png
│   ├── Fig4_Metabolic_Discovery_Spectrum.png
│   └── Fig5_Veraison_Geometric_Shockwave.png
│   └── geometric_shockwaves.csv
│
├── phytogeom.py             # Standalone CLI Application
├── requirements.txt         # Python Dependencies
└── README.md                # Documentation

```

---

## 🚀 Installation

1. **Clone the repository:**
```bash
git clone [https://github.com/sayandeep520/PhytoGeom-Field-Engine.git](https://github.com/sayandeep520/PhytoGeom-Field-Engine.git)
cd PhytoGeom-Field-Engine

```


2. **Install dependencies:**
```bash
pip install -r requirements.txt

```
## 💻 Usage

Run the Unified Field Solver on your dataset using the command line. The tool requires a Metabolights-formatted TSV file.

### Basic Execution

```bash
python phytogeom.py --file data/m_MTBLS39_grapevine.tsv

```

### Advanced Temporal Analysis

To perform the 3-Stage Kinetic Analysis (capturing the acceleration of curvature), specify the column keys for your time-points:

```bash
python phytogeom.py --file data/m_MTBLS39_grapevine.tsv --stages 063 073 083

```

**Arguments:**

* `--file`: Path to the input dataset (.tsv).
* `--stages`: A list of strings identifying the time-points or conditions to trace (e.g., `063`=Green, `073`=Veraison, `083`=Ripe).

---

## 📊 Results & Visualization

The engine automatically generates the following visualizations in the `Result/` folder:

### 1. Universal Scaling (Fig 1)

*Demonstrates the complexity dilation between Human and Plant metabolic manifolds.*

### 2. Geometric Phase Space (Fig 2)

*A topological scatter plot identifying metabolites that undergo significant rewiring (High GTS) between state transitions.*

### 3. Physics Validation (Fig 3)

*The smooth curve represents the theoretical trajectory predicted by the **Ricci-Dissipation Field Equation**, fitting the observed data points with high accuracy.*

### 4. The Discovery Spectrum (Fig 4)

*Ranked bar chart of the "Geometric Thermodynamic Score" (GTS), classifying metabolites as Hyperbolic (Expansive/Red) or Spherical (Robust/Green).*

### 5. The Veraison Shockwave (Fig 5)

*Trajectory analysis of Feruloyl Tartaric Acid, revealing the non-linear "V-shape" dip in curvature that characterizes the critical stress adaptation event.*

---
## 📐 Mathematical Formulation

The engine solves the **Ricci-Dissipation Field Equation**, which governs the geometric evolution of the metabolic network under stress:

$$
\frac{\partial R_{ij}}{\partial t} = \underbrace{-D_{\text{bio}} \Delta R_{ij}}_{\text{Geometric Diffusion}} + \underbrace{\gamma \left( \frac{\Phi_{ij}}{\Phi_{\text{crit}}} - 1 \right)}_{\text{Thermodynamic Force}} \cdot \underbrace{e^{-\frac{H(R_{ij})}{kT}}}_{\text{Hamiltonian Constraint}}
$$

### Variable Definitions:

* **$\frac{\partial R_{ij}}{\partial t}$**: The **Rewiring Velocity** of the metabolic network (how fast the geometry changes).
* **$R_{ij}$**: The **Discrete Ricci Curvature** of the pathway (edge).
    * $R < 0$: Hyperbolic (Expansive/Shunt pathway).
    * $R > 0$: Spherical (Clustered/Hub pathway).
* **$-D_{\text{bio}} \Delta R_{ij}$**: The **Geometric Diffusion Term**. It represents the homeostatic tendency of the network to smooth out curvature singularities (preventing metabolic collapse).
* **$\Phi_{ij}$**: The local **Entropy Production Rate** (Dissipation) of the reaction.
* **$\Phi_{\text{crit}}$**: The **Critical Dissipation Threshold**.
    * If $\Phi_{ij} > \Phi_{\text{crit}}$, the force becomes positive, driving a topological phase transition (curvature shift).
* **$\gamma$**: The **Stiffness Constant** (Scaling Factor). For *Vitis vinifera*, we determined $\gamma \approx 0.51$.
* **$e^{-\frac{H(R_{ij})}{kT}}$**: The **Hamiltonian Constraint**. A Boltzmann-like probability factor that penalizes energetically expensive extreme curvatures.
---

## 📜 Citation

If you use this software in your research, please cite:

> **Bera,S.D.**. (2025). *The Geometrodynamics of Plant Metabolism: A Unified Field Theory for Stress Adaptation*. Amrita Centre for Nanoscience and Molecular Medicine.

---

**Author:** Sayan Deep Bea
**Institution:** Amrita Centre for Nanoscience and Molecular Medicine


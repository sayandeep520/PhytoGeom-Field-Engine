
# PhytoGeom_Engine_v1.2.py - Comprehensive Version

# --- Imports ---
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import networkx as nx
from GraphRicciCurvature.FormanRicci import FormanRicci
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

# --- Block 1: Environment & Physics Patch ---
def setup_env():
    print("⏳ Stage 1: Installing Differential Geometry & HPC Libraries...")
    packages = ["POT", "GraphRicciCurvature", "networkx", "pandas", "numpy", "matplotlib", "seaborn", "scipy"]
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q"] + packages)

    print("⏳ Stage 2: Applying NumPy 2.0+ Compatibility Patch...")
    # Essential fix for bioinformatics legacy tools in Google Colab
    for attr in ["bool", "float", "int"]:
        if not hasattr(np, attr): setattr(np, attr, getattr(__builtins__, attr))

    print("⏳ Stage 3: Initializing 'Dark Matter' Visual Engine...")
    plt.style.use('dark_background')
    plt.rcParams['figure.facecolor'] = '#0e1117'
    plt.rcParams['axes.facecolor'] = '#12151c'
    plt.rcParams['grid.color'] = '#2a2e35'

    print("\n✅ PHYTOGEOM v1.1 ENVIRONMENT READY")

# --- Block 2: Complete PhytoGeom Methods Suite ---
class PhytoGeomEngine:
    def __init__(self, file_path, threshold=0.45):
        self.file_path = file_path
        self.threshold = threshold
        self.clean_df = None

    def ingest_data(self):
        """Robust ingestion for MetaboLights datasets."""
        print(f"📡 Ingesting: {self.file_path}")
        df = pd.read_csv(self.file_path, sep='\t')

        meta = ['database_identifier', 'chemical_formula', 'smiles', 'inchi',
                'metabolite_identification', 'mass_to_charge', 'retention_time',
                'taxid', 'species', 'database', 'reliability', 'uri', 'search_engine']

        data_cols = [c for c in df.columns if c not in meta and not c.startswith('smallmolecule_')]
        df_num = df[data_cols].copy()

        for col in df_num.columns:
            df_num[col] = pd.to_numeric(df_num[col].astype(str).str.replace(',', ''), errors='coerce')

        df_num.index = df['metabolite_identification'].fillna('Unknown')
        df_num = df_num.dropna(how='all').fillna(0)

        if not df_num.index.is_unique:
            print(f"⚠️ Merging {df_num.index.duplicated().sum()} duplicate metabolites...")
            df_num = df_num.groupby(level=0).mean()

        self.clean_df = df_num
        print(f"✅ Data Ready: {self.clean_df.shape[0]} Unique Metabolites.")
        return self.clean_df

    def build_manifold(self, target_df):
        """Constructs the Weighted Ricci Manifold for a data slice."""
        corr = target_df.T.corr(method='spearman').abs().fillna(0)

        G = nx.Graph()
        nodes = corr.columns
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                w = corr.iloc[i, j]
                if w > self.threshold:
                    G.add_edge(nodes[i], nodes[j], weight=w)

        if G.number_of_nodes() == 0:
            return {}, G

        mapping = {node: i for i, node in enumerate(G.nodes())}
        inv_mapping = {i: node for node, i in mapping.items()}
        G_int = nx.relabel_nodes(G, mapping)

        frc = FormanRicci(G_int, weight='weight')
        frc.compute_ricci_curvature()

        edge_curv = nx.get_edge_attributes(frc.G, "formanCurvature")
        node_R = {}
        for n_int in frc.G.nodes():
            incident = frc.G.edges(n_int)
            vals = [edge_curv.get((u, v), edge_curv.get((v, u), 0)) for u, v in incident]
            node_R[inv_mapping[n_int]] = np.mean(vals) if vals else 0

        return node_R, G

    def analyze_shockwave(self, key_a, key_b):
        """Solves the Geometric-Thermodynamic transition between two states."""
        cols_a = [c for c in self.clean_df.columns if key_a in c]
        cols_b = [c for c in self.clean_df.columns if key_b in c]

        if not cols_a or not cols_b:
            raise ValueError(f"Keywords '{key_a}' or '{key_b}' not found in data.")

        print(f"🕒 Analyzing Shockwave: {key_a} vs {key_b}")
        r_a, _ = self.build_manifold(self.clean_df[cols_a])
        r_b, _ = self.build_manifold(self.clean_df[cols_b])

        # Calculate Local Dissipation (Phi)
        phi_a = {n: float(self.clean_df[cols_a].loc[n].std() * np.log1p(self.clean_df[cols_a].loc[n].mean())) for n in r_a}
        phi_b = {n: float(self.clean_df[cols_b].loc[n].std() * np.log1p(self.clean_df[cols_b].loc[n].mean())) for n in r_b}

        # Manifold Alignment & Delta Calculation
        m = pd.DataFrame({'R_a': r_a, 'R_b': r_b, 'Phi_a': phi_a, 'Phi_b': phi_b}).dropna()
        m['Delta_R'] = m['R_b'] - m['R_a']
        m['Delta_Phi'] = m['Phi_b'] - m['Phi_a']
        # GTS: Geometric Tension Score (Measures non-linear structural acceleration)
        m['GTS'] = np.abs(m['Delta_R']) * (m['Phi_b'] / (m['Phi_a'] + 1e-6))

        return m.sort_values(by='GTS', ascending=False)

# --- Universal Scaling & Kingdom Comparison Function (from Block 4) ---
def solve_universal_scaling():
    print("⚖️ Initiating Universal Metric Comparison (Plant vs. Human)...")

    # 1. DEFINE DATA SOURCES
    PLANT_FILE = 'm_MTBLS39_the_plasticity_of_the_grapevine_berry_transcriptome_metabolite_profiling_mass_spectrometry_v2_maf.tsv'
    HUMAN_FILE = 'm_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'

    if not os.path.exists(PLANT_FILE) or not os.path.exists(HUMAN_FILE):
        print("❌ Error: Missing datasets. Ensure both MTBLS39 and MTBLS1 are in your sidebar.")
        return

    # --- LOCAL HELPER: Solves the 'AttributeError' by providing stats locally ---
    def get_stats_locally(curvatures, df):
        phi = {n: float(df.loc[n].std() * np.log1p(df.loc[n].mean())) for n in curvatures.keys()}
        res = pd.DataFrame({'R': curvatures, 'Phi': phi}).dropna()
        return res

    # 2. PROCESS PLANT MANIFOLD (Grapevine)
    print("🌿 Processing Plant Manifold...")
    eng_p = PhytoGeomEngine(PLANT_FILE)
    eng_p.ingest_data()
    r_p_map, _ = eng_p.build_manifold(eng_p.clean_df)
    stats_p = get_stats_locally(r_p_map, eng_p.clean_df)

    # 3. PROCESS HUMAN MANIFOLD (NMR)
    print("🧬 Processing Human Manifold...")
    eng_h = PhytoGeomEngine(HUMAN_FILE)
    eng_h.ingest_data()
    r_h_map, _ = eng_h.build_manifold(eng_h.clean_df)
    stats_h = get_stats_locally(r_h_map, eng_h.clean_df)

    # 4. COMPUTE THE UNIVERSAL SCALING FACTOR (Lambda)
    mean_r_plant = stats_p['R'].mean()
    mean_r_human = stats_h['R'].mean()
    lambda_factor = mean_r_plant / mean_r_human

    # 5. EXTREMELY HIGH-DEFINITION VISUALIZATION
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))

    # Plot A: Geometric Divergence (Density Overlap)
    sns.kdeplot(stats_p['R'], fill=True, color="#00ffcc", lw=3, label=f"Plant (Avg R: {mean_r_plant:.2f})", ax=ax[0])
    sns.kdeplot(stats_h['R'], fill=True, color="#ff0055", lw=3, label=f"Human (Avg R: {mean_r_human:.2f})", ax=ax[0])
    ax[0].set_title("A. Kingdom Manifold Divergence", fontsize=16, fontweight='bold')
    ax[0].set_xlabel("Ricci Curvature (R)")
    ax[0].legend()
    ax[0].grid(alpha=0.1)

    # Plot B: Universal Scaling Law (Log-Log)
    ax[1].scatter(stats_p['Phi'], np.abs(stats_p['R']), color="#00ffcc", alpha=0.4, s=80, label="Plant Nodes")
    ax[1].scatter(stats_h['Phi'], np.abs(stats_h['R']), color="#ff0055", alpha=0.4, s=80, label="Human Nodes")
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_title(r"B. Universal Scaling Law: Dissipation ($\Phi$) vs $|R|$ ")
    ax[1].set_xlabel(r"Dissipation ($\Phi$)")
    ax[1].set_ylabel(r"Absolute Curvature $|R|$ ")
    ax[1].grid(True, which="both", ls="-", alpha=0.05)
    ax[1].legend()

    plt.tight_layout()
    plt.show()

    print("\n" + "="*65)
    print("📐 UNIVERSAL SCALING FACTOR (LAMBDA) IDENTIFIED")
    print("="*65)
    print(f"Plant systemic curvature scale: {mean_r_plant:.4f}")
    print(f"Human systemic curvature scale: {mean_r_human:.4f}")
    print(f"CALCULATED LAMBDA (Λ): {lambda_factor:.4f}")
    print("-" * 65)
    print("💡 THEORETICAL INSIGHT:")
    if lambda_factor < 1:
        print(f"Conclusion: Human metabolism is {1/lambda_factor:.2f}x more hyperbolic than plant.")
        print("This supports the 'Complexity Dilation' theory: vertebrate manifolds")
        print("expand negatively to manage higher information flux.")
    else:
        print("Conclusion: Plant manifolds exhibit higher geometric clustering (redundancy).")
    print("="*65)

# --- Direct Execution (Main Block) ---
if __name__ == "__main__":
    setup_env()

    # --- Block 3: PhytoGeom Result Dashboard & 3D Visualization ---
    FILE = next((f for f in os.listdir('.') if f.endswith('.tsv')), None)

    if FILE:
        keys = ('063', '083') if 'MTBLS39' in FILE else ('ADG10003u', 'ADG19007u')

        pg = PhytoGeomEngine(FILE)
        pg.ingest_data()
        results = pg.analyze_shockwave(keys[0], keys[1])

        fig = plt.figure(figsize=(24, 11))

        ax1 = fig.add_subplot(131, projection='3d')
        x_val = (results['Phi_b'] / results['Phi_b'].median()) - 1
        y_val = results['R_b']**2
        z_val = results['Delta_R']

        sc = ax1.scatter(x_val, y_val, z_val, c=z_val, cmap='magma', s=110, edgecolors='w', alpha=0.8)
        ax1.set_title("A. 3D Metabolic Energy Landscape", fontsize=16, pad=20)
        ax1.set_xlabel(r"Stress $(\Phi)$")
        ax1.set_ylabel(r"Hamiltonian $H(R)$")
        ax1.set_zlabel(r"Geometric Shift $(\Delta R)$")
        plt.colorbar(sc, ax=ax1, shrink=0.5, aspect=10, label='Delta R')

        ax2 = fig.add_subplot(132)
        sns.scatterplot(data=results, x='Delta_Phi', y='Delta_R', size='GTS', hue='GTS',
                        palette='flare', sizes=(60, 800), alpha=0.6, ax=ax2)
        ax2.set_title("B. Geometric Phase Transition Map", fontsize=16)
        ax2.axhline(0, color='white', ls='--', alpha=0.3)
        ax2.axvline(0, color='white', ls='--', alpha=0.3)
        ax2.set_xlabel(r"Change in Dissipation ($\Delta \Phi$)")
        ax2.set_ylabel(r"Change in Curvature ($\Delta R$)")

        ax3 = fig.add_subplot(133)
        top_15 = results.head(15).iloc[::-1]
        ax3.barh(top_15.index, top_15['GTS'], color=sns.color_palette("viridis", 15))
        ax3.set_title("C. Top 15 Geometric Shockwave Drivers", fontsize=16)
        ax3.set_xlabel("Geometric Tension Score (GTS)")
        ax3.grid(axis='x', alpha=0.1)

        plt.tight_layout()
        plt.show()

        print("\n" + "="*60)
        print("⚛️ PHYTOGEOM v1.1 FINAL DIAGNOSTIC REPORT")
        print("="*60)
        print(f"SOURCE DATA: {FILE}")
        print(f"GEOMETRIC HUB: {results['GTS'].idxmax()}")
        print(f"SYSTEMIC REWIRING (MEAN ΔR): {results['Delta_R'].mean():.4f}")
        print("-" * 60)
        print("TOP 5 LEAPER BIOMARKERS:")
        for i, name in enumerate(results.index[:5], 1):
            print(f"  {i}. {name} (GTS: {results.loc[name, 'GTS']:.2f})")
        print("="*60)
    else:
        print("❌ ERROR: No .tsv file found. Please upload your MetaboLights file to the sidebar.")

    # --- Block 4: Universal Scaling & Kingdom Comparison Execution ---
    solve_universal_scaling()

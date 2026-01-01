
# PhytoGeom_Engine_v1.2.py

# --- Imports ---
import os
import pandas as pd
import numpy as np
import networkx as nx
from GraphRicciCurvature.FormanRicci import FormanRicci
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

# --- PhytoGeomEngine Class (from Block 2) ---
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

# --- Universal Scaling & Kingdom Comparison Function (from Block 4, modified to be a standalone function) ---
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
    # Lambda identifies the 'Curvature Ratio' between kingdoms
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

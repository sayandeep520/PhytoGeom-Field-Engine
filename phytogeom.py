"""
PHYTOGEOM: Metabolic Geometrodynamics Engine
A Unified Field Theory Solver for Plant Stress Adaptation.

Author: Sayan Deep Bera
Institution: Amrita Centre for Nanoscience and Molecular Medicine
Version: 1.0.0
"""

import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from GraphRicciCurvature.FormanRicci import FormanRicci
import os
import sys

# --- 1. THE CORE ENGINE ---
class PhytoGeomEngine:
    def __init__(self, threshold=0.25):
        self.threshold = threshold
        self.clean_df = None

    def load_data(self, file_path):
        """Ingests and cleans MTBLS datasets."""
        print(f"🔄 Loading Data from: {file_path}")
        try:
            df_raw = pd.read_csv(file_path, sep='\t')
            
            # Smart Column Detection (MN/CS/AM for Grapevine, ADG for Human)
            meta = ['database_identifier', 'chemical_formula', 'smiles', 'inchi', 'metabolite_identification']
            sample_cols = [c for c in df_raw.columns if c not in meta and not c.startswith('smallmolecule')]
            
            # Clean
            df_clean = df_raw[['metabolite_identification'] + sample_cols].copy()
            df_clean.rename(columns={'metabolite_identification': 'Metabolite'}, inplace=True)
            
            # Numeric Conversion
            for col in sample_cols:
                if df_clean[col].dtype == object:
                    df_clean[col] = df_clean[col].astype(str).str.replace(',', '').apply(pd.to_numeric, errors='coerce')
            
            # Transform
            self.clean_df = df_clean.groupby('Metabolite').mean().fillna(0)
            self.clean_df = np.log1p(self.clean_df) # Log Transform for Scaling
            
            print(f"✅ Data Loaded: {self.clean_df.shape} Matrix")
            return True
        except Exception as e:
            print(f"❌ Error Loading Data: {e}")
            return False

    def compute_curvature(self, time_point_key):
        """Calculates Ricci Curvature for a specific time point."""
        cols = [c for c in self.clean_df.columns if time_point_key in c]
        if not cols:
            print(f"⚠️ Warning: No samples found for key '{time_point_key}'")
            return {}, None

        sub_df = self.clean_df[cols]
        corr = sub_df.T.corr(method='spearman').abs().fillna(0)
        
        G = nx.Graph()
        nodes = corr.columns
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if corr.iloc[i, j] > self.threshold:
                    G.add_edge(nodes[i], nodes[j], weight=corr.iloc[i, j])

        # Integer Mapping Patch for GraphRicciCurvature
        mapping = {node: i for i, node in enumerate(G.nodes())}
        inv_mapping = {i: node for node, i in mapping.items()}
        G_int = nx.relabel_nodes(G, mapping)
        
        if len(G_int.edges()) > 0:
            frc = FormanRicci(G_int, weight='weight')
            frc.compute_ricci_curvature()
            edge_curv = nx.get_edge_attributes(frc.G, "formanCurvature")
            
            node_R = {}
            for n_int in frc.G.nodes():
                incident = frc.G.edges(n_int)
                vals = [edge_curv.get((u, v), edge_curv.get((v, u), 0)) for u, v in incident]
                node_R[inv_mapping[n_int]] = np.mean(vals) if vals else 0
            return node_R, G
        return {}, G

# --- 2. THE PHYSICS SOLVER ---
def ricci_flow_equation(R, t, force, gamma=0.5, D=0.1):
    """The Ricci-Dissipation Field Equation."""
    geometric_hamiltonian = R**2
    boltzmann = np.exp(-geometric_hamiltonian)
    dR_dt = -D * R + gamma * (force - 1) * boltzmann
    return dR_dt

def simulate_trajectory(metabolite, path):
    """Simulates the physics of the adaptation."""
    t = np.linspace(0, 100, 100)
    # Estimate force required to bridge the gap
    force = 1.5 if path[-1] > path[0] else 0.5
    solution = odeint(ricci_flow_equation, path[0], t, args=(force,))
    return t, solution.flatten()

# --- 3. MAIN EXECUTION ---
def main():
    parser = argparse.ArgumentParser(description="Metabolic Geometrodynamics Solver")
    parser.add_argument('--file', type=str, required=True, help='Path to the TSV dataset')
    parser.add_argument('--stages', nargs='+', default=['063', '073', '083'], help='List of time-point keys (e.g., 063 073 083)')
    args = parser.parse_args()

    engine = PhytoGeomEngine(threshold=0.25)
    if not engine.load_data(args.file):
        return

    print("\n🌊 INITIATING 3-STAGE KINETIC ANALYSIS")
    curvature_data = {}
    
    # 1. Compute Curvature
    for stage in args.stages:
        print(f"   📍 Mapping Manifold for Stage: {stage}...")
        R, _ = engine.compute_curvature(stage)
        curvature_data[stage] = R

    # 2. Find Accelerators
    common = set.intersection(*[set(d.keys()) for d in curvature_data.values()])
    print(f"   🔗 Analyzing {len(common)} persistent pathways...")

    results = []
    for m in common:
        path = [curvature_data[s][m] for s in args.stages]
        # Acceleration ~ |R_final - 2*R_mid + R_start|
        accel = abs(path[2] - 2*path[1] + path[0])
        results.append({'Metabolite': m, 'R_start': path[0], 'R_mid': path[1], 'R_final': path[2], 'Acceleration': accel})

    df = pd.DataFrame(results).sort_values(by='Acceleration', ascending=False)
    
    # 3. Report
    print("\n🏆 TOP 5 METABOLIC DRIVERS (Geometric Shockwaves)")
    print(df.head(5).to_string(index=False))
    
    # 4. Save
    if not os.path.exists('output'): os.makedirs('output')
    df.to_csv('output/geometric_shockwaves.csv', index=False)
    print("\n✅ Results saved to output/geometric_shockwaves.csv")

    # 5. Plot Top Driver
    top_driver = df.iloc[0]
    name = top_driver['Metabolite']
    path = [top_driver['R_start'], top_driver['R_mid'], top_driver['R_final']]
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(path)), path, marker='o', color='#00ff9d', linewidth=2, label='Observed Path')
    
    # Add Physics Simulation Overlay
    t_sim, r_sim = simulate_trajectory(name, path)
    # Scale time to match plot x-axis (0 to 2)
    t_scaled = t_sim / 100 * 2 
    plt.plot(t_scaled, r_sim, '--', color='white', alpha=0.5, label='Field Eq. Theory')
    
    plt.title(f"Geometric Evolution of {name}")
    plt.ylabel("Ricci Curvature")
    plt.xticks(range(len(args.stages)), args.stages)
    plt.legend()
    plt.savefig(f'output/trajectory_{name}.png')
    print(f"✅ Trajectory plot saved to output/trajectory_{name}.png")

if __name__ == "__main__":
    main()

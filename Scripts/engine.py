import pandas as pd
import numpy as np
import networkx as nx
from GraphRicciCurvature.FormanRicci import FormanRicci

class PhytoGeomEngine:
    def __init__(self, file_path=None, threshold=0.45):
        self.file_path = file_path
        self.threshold = threshold
        self.clean_df = None

    def ingest_data(self, df=None):
        if df is None:
            df = pd.read_csv(self.file_path, sep='\t')

        # Metadata filter
        meta = ['database_identifier', 'chemical_formula', 'smiles', 'inchi',
                'metabolite_identification', 'mass_to_charge', 'retention_time']
        data_cols = [c for c in df.columns if c not in meta and not c.startswith('smallmolecule_')]
        df_num = df[data_cols].copy()

        for col in df_num.columns:
            df_num[col] = pd.to_numeric(df_num[col].astype(str).str.replace(',', ''), errors='coerce')

        df_num.index = df['metabolite_identification'].fillna('Unknown')
        self.clean_df = df_num.dropna(how='all').fillna(0).groupby(level=0).mean()
        print(f"✅ Data Ready: {self.clean_df.shape[0]} Unique Metabolites.")
        return self.clean_df

    def build_manifold(self, target_df):
        corr = target_df.T.corr(method='spearman').abs().fillna(0)
        G = nx.Graph()
        nodes = corr.columns
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                w = corr.iloc[i, j]
                if w > self.threshold:
                    G.add_edge(nodes[i], nodes[j], weight=w)

        mapping = {node: i for i, node in enumerate(G.nodes())}
        inv_mapping = {i: node for node, i in mapping.items()}
        G_int = nx.relabel_nodes(G, mapping)

        frc = FormanRicci(G_int, weight='weight')
        frc.compute_ricci_curvature()

        node_R = {}
        edge_curv = nx.get_edge_attributes(frc.G, "formanCurvature")
        for n_int in frc.G.nodes():
            incident = frc.G.edges(n_int)
            vals = [edge_curv.get((u, v), edge_curv.get((v, u), 0)) for u, v in incident]
            node_R[inv_mapping[n_int]] = np.mean(vals) if vals else 0
        return node_R, G

    def analyze_shockwave(self, key_a, key_b):
        cols_a = [c for c in self.clean_df.columns if key_a in c]
        cols_b = [c for c in self.clean_df.columns if key_b in c]
        r_a, _ = self.build_manifold(self.clean_df[cols_a])
        r_b, _ = self.build_manifold(self.clean_df[cols_b])

        phi_a = {n: float(self.clean_df[cols_a].loc[n].std() * np.log1p(self.clean_df[cols_a].loc[n].mean())) for n in r_a}
        phi_b = {n: float(self.clean_df[cols_b].loc[n].std() * np.log1p(self.clean_df[cols_b].loc[n].mean())) for n in r_b}

        m = pd.DataFrame({'R_a': r_a, 'R_b': r_b, 'Phi_a': phi_a, 'Phi_b': phi_b}).dropna()
        m['GTS'] = np.abs(m['R_b'] - m['R_a']) * (m['Phi_b'] / (m['Phi_a'] + 1e-6))
        return m.sort_values(by='GTS', ascending=False)

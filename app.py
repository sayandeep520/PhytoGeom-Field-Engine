import streamlit as st
import pandas as pd
import numpy as np
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
from scipy.integrate import odeint
from GraphRicciCurvature.FormanRicci import FormanRicci

# --- CONFIGURATION ---
st.set_page_config(page_title="PhytoGeom Web Suite", layout="wide", page_icon="🌿")

# --- CSS STYLING ---
st.markdown("""
<style>
    .reportview-container { background: #0e1117; }
    h1 { color: #00ff9d; }
    h2, h3 { color: white; }
    .stButton>button { background-color: #00ff9d; color: black; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

# --- 1. THE ENGINE (Cached for Speed) ---
@st.cache_data
def load_data(uploaded_file):
    try:
        df_raw = pd.read_csv(uploaded_file, sep='\t')
        meta = ['database_identifier', 'chemical_formula', 'smiles', 'inchi', 'metabolite_identification']
        sample_cols = [c for c in df_raw.columns if c not in meta and not c.startswith('smallmolecule')]
        df_clean = df_raw[['metabolite_identification'] + sample_cols].copy()
        df_clean.rename(columns={'metabolite_identification': 'Metabolite'}, inplace=True)
        for col in sample_cols:
            if df_clean[col].dtype == object:
                df_clean[col] = df_clean[col].astype(str).str.replace(',', '').apply(pd.to_numeric, errors='coerce')
        df_final = df_clean.groupby('Metabolite').mean().fillna(0)
        df_final = np.log1p(df_final)
        return df_final
    except Exception as e:
        return None

def compute_manifold(df, threshold):
    corr = df.T.corr(method='spearman').abs().fillna(0)
    G = nx.Graph()
    nodes = corr.columns
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if corr.iloc[i, j] > threshold:
                G.add_edge(nodes[i], nodes[j], weight=corr.iloc[i, j])
    
    mapping = {node: i for i, node in enumerate(G.nodes())}
    inv_mapping = {i: node for node, i in mapping.items()}
    G_int = nx.relabel_nodes(G, mapping)
    
    node_R = {}
    if len(G_int.edges()) > 0:
        frc = FormanRicci(G_int, weight='weight')
        frc.compute_ricci_curvature()
        edge_curv = nx.get_edge_attributes(frc.G, "formanCurvature")
        for n_int in frc.G.nodes():
            incident = frc.G.edges(n_int)
            vals = [edge_curv.get((u, v), edge_curv.get((v, u), 0)) for u, v in incident]
            node_R[inv_mapping[n_int]] = np.mean(vals) if vals else 0
            
    return node_R, G

# --- 2. THE UI ---
st.title("🌿 PhytoGeom: Unified Field Solver")
st.markdown("### The Geometrodynamics of Plant Metabolism")

# Sidebar
st.sidebar.header("1. Upload Data")
uploaded_file = st.sidebar.file_uploader("Choose a Metabolights TSV", type="tsv")

st.sidebar.header("2. Model Parameters")
threshold = st.sidebar.slider("Connectivity Threshold", 0.1, 0.9, 0.25)
gamma = st.sidebar.slider("Stiffness Constant (γ)", 0.1, 1.0, 0.51)
st.sidebar.markdown("---")
st.sidebar.info("Developed by Amrita Centre for Nanoscience")

if uploaded_file is not None:
    df = load_data(uploaded_file)
    if df is not None:
        st.success(f"Dataset Loaded: {df.shape[0]} Metabolites x {df.shape[1]} Samples")
        
        # Define Stages
        col1, col2 = st.columns(2)
        with col1:
            stage_a = st.text_input("Initial Stage Key (e.g., 063)", "063")
        with col2:
            stage_b = st.text_input("Final Stage Key (e.g., 083)", "083")

        if st.button("🚀 RUN UNIFIED FIELD SIMULATION"):
            with st.spinner("Solving Ricci Flow Equations..."):
                # Data Slicing
                cols_a = [c for c in df.columns if stage_a in c]
                cols_b = [c for c in df.columns if stage_b in c]
                
                if cols_a and cols_b:
                    # Computation
                    r_a, G_a = compute_manifold(df[cols_a], threshold)
                    r_b, G_b = compute_manifold(df[cols_b], threshold)
                    
                    # Analysis
                    common = set(r_a.keys()) & set(r_b.keys())
                    results = []
                    for m in common:
                        results.append({
                            'Metabolite': m,
                            'R_green': r_a[m],
                            'R_ripe': r_b[m],
                            'GTS': abs(r_b[m] - r_a[m])
                        })
                    
                    res_df = pd.DataFrame(results).sort_values(by='GTS', ascending=False)
                    
                    # --- RESULTS DASHBOARD ---
                    tab1, tab2, tab3 = st.tabs(["📉 Geometric Shockwaves", "🕸 Network Topology", "🧬 Physics Simulation"])
                    
                    with tab1:
                        st.subheader("Topological Phase Transition")
                        fig = px.scatter(res_df, x="R_green", y="R_ripe", 
                                       size="GTS", color="GTS", hover_name="Metabolite",
                                       color_continuous_scale="Turbo", size_max=40,
                                       title=f"Phase Space: {stage_a} → {stage_b}")
                        fig.add_shape(type="line", x0=-50, y0=-50, x1=50, y1=50, 
                                    line=dict(color="White", dash="dash"))
                        st.plotly_chart(fig, use_container_width=True)
                        
                        st.subheader("Top Discovery Leads")
                        st.dataframe(res_df.head(10))

                    with tab2:
                        st.subheader("Metabolic Manifold Structure")
                        # Simple metric display
                        k_avg = np.mean(list(r_b.values()))
                        st.metric("Average System Curvature (R)", f"{k_avg:.4f}")
                        st.write(f"Nodes: {len(G_b.nodes())} | Edges: {len(G_b.edges())}")
                        st.warning("3D Interactive Graph is disabled in Lite Mode to save memory.")

                    with tab3:
                        st.subheader("Ricci-Dissipation Trajectory")
                        selected_met = st.selectbox("Select Metabolite", res_df['Metabolite'].head(20))
                        
                        row = res_df[res_df['Metabolite'] == selected_met].iloc[0]
                        start, end = row['R_green'], row['R_ripe']
                        
                        # Simulate
                        t = np.linspace(0, 10, 100)
                        def model(y, t): return -0.1 * y + gamma * (1.5 if end > start else 0.5 - 1) * np.exp(-y**2)
                        y = odeint(model, start, t)
                        
                        fig_phys = px.line(x=t, y=y.flatten(), labels={'x': 'Time', 'y': 'Curvature (R)'})
                        fig_phys.add_scatter(x=[0, 10], y=[start, end], mode='markers', marker=dict(size=15, color='red'), name='Observed')
                        fig_phys.update_layout(title=f"Evolution of {selected_met}")
                        st.plotly_chart(fig_phys, use_container_width=True)

                else:
                    st.error("Could not find columns matching those keys.")

import gudhi
import matplotlib.pyplot as plt

def run_topology_audit(G, results_df):
    print("🕸️ Analyzing Manifold Persistence (Betti Numbers)...")

    # Integer mapping for GUDHI C++ engine
    node_to_int = {node: i for i, node in enumerate(G.nodes())}
    st = gudhi.SimplexTree()

    for u, v, d in G.edges(data=True):
        # Distance = 1 - Correlation (Filtration)
        st.insert([node_to_int[u], node_to_int[v]], filtration=1.0 - d.get('weight', 0.5))

    st.compute_persistence()
    betti = st.betti_numbers()

    # Visualization: Persistence Barcode
    plt.figure(figsize=(10, 3))
    gudhi.plot_persistence_barcode(st.persistence())
    plt.title("Structural Feedback Stability (Persistence Barcode)", fontsize=12)
    plt.xlabel("Filtration Value (1 - Correlation)")
    plt.tight_layout()
    plt.show()

    print(f"📊 Audit: {betti[0]} Integrated Clusters | {betti[1] if len(betti)>1 else 0} Regulatory Loops.")

    # Final Interpretation & Annotation
    # We map the physics leads to their probable biological roles
    role_map = {
        'Proline': 'Osmolyte / Structural Shield',
        'GABA': 'Stress Signaling',
        'Malate': 'TCA Intermediate / Carbon Flux',
        'Sucrose': 'Energy Resource Management',
        'b-hydroxy-n-butyrate': 'Alternative Energy Flux',
        '4-aminohippurate': 'Secondary Defense Hub'
    }

    report = results_df.copy()
    report['Biological_Role'] = report.index.map(lambda x: role_map.get(str(x), 'Discovery Lead'))

    print("\n" + "="*75)
    print("🚀 GEOMETRODYNAMIC UNIFIED FIELD REPORT: FINAL SCIENTIFIC DISCOVERY")
    print("="*75)
    cols = [c for c in ['GTS', 'R_b', 'Phi_b', 'Biological_Role'] if c in report.columns]
    print(report[cols].sort_values('GTS', ascending=False).head(10))
    print("="*75)

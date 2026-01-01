import numpy as np
import networkx as nx
from scipy.integrate import odeint

class FieldSolver:
    def __init__(self, G, node_data, results_df):
        self.G = G
        self.nodes = list(G.nodes())
        self.node_idx = {n: i for i, n in enumerate(self.nodes)}
        # Laplacian Matrix: Represents the 'network fabric'
        self.L = nx.laplacian_matrix(G, weight='weight').toarray()

        # Phi: Thermodynamic Dissipation proxy
        self.Phi = np.array([node_data[n]['std'] * np.log1p(abs(node_data[n]['mean'])) for n in self.nodes])

        # Initial State: R_b (The curvature at the height of the shockwave)
        self.R0 = np.array([results_df.loc[n, 'R_b'] if n in results_df.index else 0 for n in self.nodes])

    def field_equation(self, R, t):
        # ∂R/∂t = -D_bio * ΔR + γ * (Φ/Φ_crit - 1) * exp(-H(R)/kT)
        D_bio = 0.12   # Diffusion constant (Structural healing rate)
        gamma = 0.55   # Coupling constant (Stiffness)
        phi_crit = 1.2 # Critical threshold

        diffusion = -D_bio * (self.L @ R)
        source = gamma * (self.Phi / phi_crit - 1)
        damping = np.exp(-(0.015 * R**2)) # Hamiltonian penalty

        return diffusion + (source * damping)

    def solve(self, tau_max=10):
        print(f"🧬 Simulating Field Evolution for {len(self.nodes)} metabolites...")
        t = np.linspace(0, tau_max, 100)
        R_path = odeint(self.field_equation, self.R0, t)
        return t, R_path

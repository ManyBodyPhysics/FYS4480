import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.optim as optim


# =============================================================================
# 1. Fermionic operators in Fock space
# =============================================================================

def build_fermionic_ops(n_orb):
    """
    Build creation and annihilation operators for n_orb spin orbitals
    using Jordan–Wigner representation directly in Fock space.

    Basis: computational basis of dimension dim = 2^n_orb, with bit k
    indicating occupation of orbital k.

    Returns:
        a_dag: list of length n_orb with creation matrices (dim x dim)
        a    : list of length n_orb with annihilation matrices (dim x dim)
    """
    dim = 2 ** n_orb
    a_dag = []
    a = []

    for k in range(n_orb):
        create = np.zeros((dim, dim), dtype=np.complex128)
        annih = np.zeros((dim, dim), dtype=np.complex128)

        for b in range(dim):
            occ = (b >> k) & 1
            # phase = (-1)^(number of occupied orbitals below k)
            below = b & ((1 << k) - 1)
            phase = -1 if (below.bit_count() % 2 == 1) else 1

            # creation: if orbital k is empty
            if occ == 0:
                b_new = b | (1 << k)
                create[b_new, b] = phase

            # annihilation: if orbital k is occupied
            if occ == 1:
                b_new = b & ~(1 << k)
                annih[b_new, b] = phase

        a_dag.append(create)
        a.append(annih)

    return a_dag, a


def build_pair_ops(n_levels, a_dag, a):
    """
    Build pair creation/annihilation operators for each spatial level p.

    We assume:
        orbitals 2p = p↑, 2p+1 = p↓.

    Returns:
        P_dag: list of length n_levels with P_p^† = c_{p↑}^† c_{p↓}^†
        P    : list of length n_levels with P_p   = c_{p↓} c_{p↑}
    """
    P_dag = []
    P = []
    for p in range(n_levels):
        up = 2 * p
        dn = 2 * p + 1

        # P_p^† = c_up^† c_dn^†
        Pp_dag = a_dag[up] @ a_dag[dn]
        # P_p   = c_dn c_up
        Pp = a[dn] @ a[up]

        P_dag.append(Pp_dag)
        P.append(Pp)
    return P_dag, P


def build_hf_state(n_levels, n_pairs, dim):
    """
    Hartree–Fock reference state: first n_pairs spatial levels doubly occupied,
    others empty.

    orbitals: 2p = p↑, 2p+1 = p↓
    HF bitstring: 1 for orbitals in occupied levels, else 0.
    """
    occ_orbitals = []
    for p in range(n_pairs):
        occ_orbitals.append(2 * p)
        occ_orbitals.append(2 * p + 1)

    index = 0
    for k in occ_orbitals:
        index |= (1 << k)

    psi = np.zeros(dim, dtype=np.complex128)
    psi[index] = 1.0
    return psi


# =============================================================================
# 2. Pairing Hamiltonian and UCC-like generators
# =============================================================================

@dataclass
class PairingHamiltonian:
    """
    Pairing Hamiltonian with L spatial levels and 2L spin orbitals.

    H = sum_p eps_p (n_{p↑} + n_{p↓}) - g sum_{p,q} P_p^† P_q
    eps_p: list of length n_levels
    g: pairing strength
    """
    eps_list: list
    g: float

    def build_hamiltonian(self, a_dag, a, P_dag, P):
        n_levels = len(self.eps_list)
        n_orb = 2 * n_levels
        dim = 2 ** n_orb
        H = np.zeros((dim, dim), dtype=np.complex128)

        # Single-particle part
        for p, eps in enumerate(self.eps_list):
            for spin in range(2):
                k = 2 * p + spin
                n_op = a_dag[k] @ a[k]
                H += eps * n_op

        # Pairing part: -g sum_{p,q} P_p^† P_q
        H_pair = np.zeros_like(H)
        for p in range(n_levels):
            for q in range(n_levels):
                H_pair += P_dag[p] @ P[q]

        H += -self.g * H_pair

        # Enforce Hermiticity numerically
        H = 0.5 * (H + H.conj().T)
        return H

    @staticmethod
    def exact_ground_energy(H):
        evals, _ = np.linalg.eigh(H)
        return float(np.min(evals))


def build_ucc_pair_generators(n_levels, n_pairs, P_dag, P):
    """
    Build UCC-like pair excitation generators A_j.

    Occ levels: p = 0,...,n_pairs-1
    Virt levels: q = n_pairs,...,n_levels-1

    For each (p, q), define:
        A_{p,q} = i (P_q^† P_p - P_p^† P_q)
    which is Hermitian (A = A^†).

    Returns:
        A_list: list of Hermitian matrices (dim x dim)
    """
    n_orb = 2 * n_levels
    dim = 2 ** n_orb

    occ = list(range(n_pairs))
    virt = list(range(n_pairs, n_levels))

    A_list = []
    for p in occ:
        for q in virt:
            A = 1j * (P_dag[q] @ P[p] - P_dag[p] @ P[q])
            A = 0.5 * (A + A.conj().T)  # enforce Hermitian
            A_list.append(A)

    if len(A_list) == 0:
        # trivial case: no generators
        A_list.append(np.zeros((dim, dim), dtype=np.complex128))

    return A_list


def precompute_eigendecomp(A_list):
    """
    Precompute eigenvalues/eigenvectors of each Hermitian generator A_j.

    Returns:
        spec_list: list of (eigvals, eigvecs)
    """
    spec_list = []
    for A in A_list:
        vals, vecs = np.linalg.eigh(A)
        spec_list.append((vals, vecs))
    return spec_list


def apply_unitary_from_spec(psi, eigvals, eigvecs, dt):
    """
    Apply U_j(dt) = exp(-i dt A_j) to state psi, given eigen-decomposition
    A_j v = eigvals v.

    psi: (dim,) complex
    eigvals: (dim,)
    eigvecs: (dim, dim)
    """
    phase = np.exp(-1j * dt * eigvals)
    tmp = eigvecs.conj().T @ psi
    tmp = phase * tmp
    psi_new = eigvecs @ tmp
    return psi_new


def trotterized_ucc_state(psi0, A_spec_list, theta, n_steps):
    """
    Compute the Trotterized UCC state:

        U(θ) ≈ [ ∏_j exp(-i (θ / n_steps) A_j) ]^{n_steps} |psi0>

    psi0: initial state (dim,)
    A_spec_list: list of (eigvals, eigvecs) for each A_j
    theta: UCC amplitude
    n_steps: Trotter steps

    Returns:
        psi: final state (dim,)
    """
    if n_steps <= 0:
        return psi0.copy()

    dt = theta / float(n_steps)
    psi = psi0.copy()

    for _ in range(n_steps):
        for eigvals, eigvecs in A_spec_list:
            psi = apply_unitary_from_spec(psi, eigvals, eigvecs, dt)

    norm = np.linalg.norm(psi)
    if norm > 1e-14:
        psi /= norm
    return psi


# =============================================================================
# 3. Dataset generation: E(dt, g, theta) from Trotterized UCC
# =============================================================================

def generate_dataset_pairing_ucc(
    n_levels=3,
    n_pairs=2,
    g_values=(0.4, 0.8, 1.2, 1.6),
    theta_values=(0.2, 0.4, 0.6, 0.8, 1.0),
    n_steps_values=(1, 2, 3, 4, 6, 8),
    eps0=0.0,
    delta_eps=1.0
):
    """
    Generate dataset for pairing Hamiltonian with UCC-like Trotterized ansatz.

    System:
        L = n_levels spatial levels, 2L spin orbitals, 2L qubits.
        HF reference has n_pairs doubly occupied levels.

    For each g in g_values:
        - Build H(g) using eps_p = eps0 + p * delta_eps.
        - Compute exact ground energy E_exact(g).

        For each theta in theta_values:
            For each n_steps in n_steps_values:
                - Construct Trotterized UCC state |psi(theta, n_steps)>,
                  using the pair-excitation generators A_j.
                - Compute energy E_trotter = <psi|H|psi>.
                - Define dt = theta / n_steps.
                - Store feature x = [dt, g, theta], target y = E_trotter.

    Returns:
        X: (N_samples, 3) feature matrix
        y: (N_samples,)   Trotterized energies
        y_exact: (N_samples,) exact ground energies (per sample's g)
    """
    n_orb = 2 * n_levels
    dim = 2 ** n_orb

    # Build operators independent of g
    a_dag, a = build_fermionic_ops(n_orb)
    P_dag, P = build_pair_ops(n_levels, a_dag, a)
    A_list = build_ucc_pair_generators(n_levels, n_pairs, P_dag, P)
    A_spec_list = precompute_eigendecomp(A_list)
    psi0 = build_hf_state(n_levels, n_pairs, dim)

    # Single-particle energies
    eps_list = [eps0 + p * delta_eps for p in range(n_levels)]

    X_list = []
    y_list = []
    y_exact_list = []

    for g in g_values:
        H_model = PairingHamiltonian(eps_list=eps_list, g=g)
        H = H_model.build_hamiltonian(a_dag, a, P_dag, P)
        E_exact = PairingHamiltonian.exact_ground_energy(H)

        for theta in theta_values:
            for n_steps in n_steps_values:
                if n_steps <= 0:
                    continue
                dt = theta / float(n_steps)

                psi = trotterized_ucc_state(psi0, A_spec_list, theta, n_steps)
                E_trotter = np.real(psi.conj().T @ (H @ psi))

                X_list.append([dt, g, theta])
                y_list.append(E_trotter)
                y_exact_list.append(E_exact)

    X = np.array(X_list, dtype=np.float32)
    y = np.array(y_list, dtype=np.float32)
    y_exact = np.array(y_exact_list, dtype=np.float32)
    return X, y, y_exact


# =============================================================================
# 4. Regressors: MLP and PMM
# =============================================================================

class MLPEnergyRegressor(nn.Module):
    """
    Simple MLP regressor for E(dt, g, theta).
    """
    def __init__(self, input_dim=3, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


class PMMEnergyRegressor(nn.Module):
    """
    Parametric Matrix Model (PMM) regressor for energy:

        x in R^p (p=3: dt, g, theta)

    Primary matrix:
        H(x) = H0 + sum_l x_l H_l   (complex Hermitian)
    Secondary matrix S_E (complex Hermitian).

    Energy:
        E(x) = bias + sum_{i=1}^r Re( v_i^† S_E v_i ),
    where v_i are eigenvectors corresponding to the lowest r eigenvalues
    of H(x).
    """
    def __init__(self, input_dim: int, latent_dim: int = 8, r: int = 3):
        super().__init__()
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.r = r

        # Primary matrices (complex)
        self.H0 = nn.Parameter(
            0.1 * (torch.randn(latent_dim, latent_dim, dtype=torch.complex64))
        )
        self.H_inputs = nn.Parameter(
            0.1 * (torch.randn(input_dim, latent_dim, latent_dim, dtype=torch.complex64))
        )

        # Secondary matrix
        self.S_E = nn.Parameter(
            0.1 * (torch.randn(latent_dim, latent_dim, dtype=torch.complex64))
        )

        self.bias = nn.Parameter(torch.zeros(1, dtype=torch.float32))

    @staticmethod
    def hermitian(A: torch.Tensor) -> torch.Tensor:
        return 0.5 * (A + A.conj().transpose(-1, -2))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (batch, input_dim), real tensor
        returns:
            E: (batch,), real tensor
        """
        batch_size = x.shape[0]
        x_real = x.to(torch.float32)

        H0_herm = self.hermitian(self.H0)
        H_inputs_herm = self.hermitian(self.H_inputs)
        S_E_herm = self.hermitian(self.S_E)

        # Build H(x) per sample
        H = H0_herm.unsqueeze(0).expand(batch_size, -1, -1).clone().to(torch.complex64)
        for l in range(self.input_dim):
            coef = x_real[:, l].view(-1, 1, 1).to(torch.complex64)
            H = H + coef * H_inputs_herm[l]
        H = self.hermitian(H)

        # Eigen-decomposition
        eigvals, eigvecs = torch.linalg.eigh(H)  # (b, n), (b, n, n)
        v = eigvecs[:, :, :self.r]               # (b, n, r)

        # Energy
        E_complex = torch.zeros(batch_size, dtype=torch.complex64, device=x.device)
        for i in range(self.r):
            v_i = v[:, :, i]  # (b, n)
            tmp = torch.einsum("bi,ij,bj->b", v_i.conj(), S_E_herm, v_i)
            E_complex += tmp

        E = E_complex.real + self.bias
        return E


def train_regressor(
    model,
    X: np.ndarray,
    y: np.ndarray,
    lr: float = 1e-3,
    n_epochs: int = 400,
    batch_size: int = 64,
    name: str = "model"
):
    """
    Generic training loop for both MLP and PMM regressors.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    X_t = torch.tensor(X, dtype=torch.float32, device=device)
    y_t = torch.tensor(y, dtype=torch.float32, device=device)

    N = X.shape[0]
    perm = np.random.permutation(N)
    N_train = int(0.8 * N)
    train_idx = perm[:N_train]
    val_idx = perm[N_train:]

    X_train = X_t[train_idx]
    y_train = y_t[train_idx]
    X_val = X_t[val_idx]
    y_val = y_t[val_idx]

    model = model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()

    print(f"\nTraining {name}...")
    for epoch in range(1, n_epochs + 1):
        model.train()
        perm_torch = torch.randperm(N_train, device=device)
        total_loss = 0.0
        n_batches = 0

        for i in range(0, N_train, batch_size):
            idx = perm_torch[i:i+batch_size]
            xb = X_train[idx]
            yb = y_train[idx]

            optimizer.zero_grad()
            yp = model(xb)
            loss = loss_fn(yp, yb)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        train_loss = total_loss / max(1, n_batches)
        model.eval()
        with torch.no_grad():
            yv_pred = model(X_val)
            val_loss = loss_fn(yv_pred, y_val).item()

        if epoch == 1 or epoch % 50 == 0:
            print(f"[{name}] Epoch {epoch:4d}: "
                  f"train MSE = {train_loss:.4e}, val MSE = {val_loss:.4e}")

    return model


# =============================================================================
# 5. Main: full pipeline + extrapolation plot
# =============================================================================

def main():
    # Parameters for the pairing model and dataset
    n_levels = 3   # spatial levels -> 2 * n_levels spin orbitals
    n_pairs = 2    # number of occupied pairs in HF reference

    g_values = (0.4, 0.8, 1.2, 1.6)
    theta_values = (0.2, 0.4, 0.6, 0.8, 1.0)
    n_steps_values = (1, 2, 3, 4, 6, 8)

    print("Generating dataset from Trotterized UCC pair ansatz...")
    X, y, y_exact = generate_dataset_pairing_ucc(
        n_levels=n_levels,
        n_pairs=n_pairs,
        g_values=g_values,
        theta_values=theta_values,
        n_steps_values=n_steps_values,
        eps0=0.0,
        delta_eps=1.0
    )

    print(f"Dataset size: {X.shape[0]} samples")
    print(f"Feature dimension: {X.shape[1]} (dt, g, theta)")

    # Train MLP regressor
    mlp = MLPEnergyRegressor(input_dim=X.shape[1], hidden_dim=64)
    mlp = train_regressor(mlp, X, y, lr=1e-3, n_epochs=400,
                          batch_size=64, name="MLP")

    # Train PMM regressor
    pmm = PMMEnergyRegressor(input_dim=X.shape[1], latent_dim=8, r=3)
    pmm = train_regressor(pmm, X, y, lr=1e-3, n_epochs=400,
                          batch_size=64, name="PMM")

    # -------------------------------------------------------------------------
    # Extrapolation test for a chosen (g*, theta*)
    # -------------------------------------------------------------------------
    print("\nPerforming zero-Trotter-step extrapolation test...")

    # Choose test parameters not necessarily in training grid
    g_star = 1.0
    theta_star = 0.7

    # Rebuild operators and Hamiltonian for this g*
    n_orb = 2 * n_levels
    dim = 2 ** n_orb
    a_dag, a = build_fermionic_ops(n_orb)
    P_dag, P = build_pair_ops(n_levels, a_dag, a)
    A_list = build_ucc_pair_generators(n_levels, n_pairs, P_dag, P)
    A_spec_list = precompute_eigendecomp(A_list)
    psi0 = build_hf_state(n_levels, n_pairs, dim)
    eps_list = [0.0 + p * 1.0 for p in range(n_levels)]
    H_model_star = PairingHamiltonian(eps_list=eps_list, g=g_star)
    H_star = H_model_star.build_hamiltonian(a_dag, a, P_dag, P)
    E_exact_star = PairingHamiltonian.exact_ground_energy(H_star)

    # Evaluate Trotterized energies for multiple n_steps
    n_steps_eval = np.array([1, 2, 3, 4, 6, 8, 10, 12, 16], dtype=int)
    dt_eval = theta_star / n_steps_eval.astype(np.float32)
    E_trot = []
    for n_steps in n_steps_eval:
        psi = trotterized_ucc_state(psi0, A_spec_list, theta_star, n_steps)
        E = np.real(psi.conj().T @ (H_star @ psi))
        E_trot.append(E)
    E_trot = np.array(E_trot, dtype=np.float32)

    # Predict energies via MLP and PMM
    X_eval = np.stack([dt_eval,
                       g_star * np.ones_like(dt_eval),
                       theta_star * np.ones_like(dt_eval)], axis=1)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    mlp.eval()
    pmm.eval()
    with torch.no_grad():
        X_eval_t = torch.tensor(X_eval, dtype=torch.float32, device=device)
        E_mlp = mlp(X_eval_t).cpu().numpy()
        E_pmm = pmm(X_eval_t).cpu().numpy()

    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(8, 5))
    plt.scatter(dt_eval, E_trot, label="Trotterized UCC energy", marker="o")
    plt.plot(dt_eval, E_mlp, label="MLP prediction", linestyle="--")
    plt.plot(dt_eval, E_pmm, label="PMM prediction", linestyle="-.")

    plt.axhline(E_exact_star, color="k", linestyle=":", label="Exact ground energy")

    plt.xlabel(r"Trotter step size $dt = \theta / n_{\mathrm{steps}}$")
    plt.ylabel("Energy")
    plt.title(f"Zero-Trotter-step extrapolation (g={g_star:.2f}, theta={theta_star:.2f})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Print some numerical comparisons for smallest dt
    print("\nSample numerical comparisons (smallest dt values):")
    order = np.argsort(dt_eval)
    for idx in order[:5]:
        print(f"n_steps = {n_steps_eval[idx]:2d}, "
              f"dt = {dt_eval[idx]:.4f} | "
              f"E_trot = {E_trot[idx]: .6f}, "
              f"E_MLP = {E_mlp[idx]: .6f}, "
              f"E_PMM = {E_pmm[idx]: .6f}, "
              f"E_exact = {E_exact_star: .6f}")


if __name__ == "__main__":
    main()

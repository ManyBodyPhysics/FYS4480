import numpy as np
class HartreeFock:
    def __init__(self, num_electrons, num_orbitals):
        self.num_electrons = num_electrons
        self.num_orbitals = num_orbitals
        self.h = np.random.rand(num_orbitals, num_orbitals)  # One-electron integrals
        self.coulomb = np.random.rand(num_orbitals, num_orbitals, num_orbitals, num_orbitals)  # Two-electron integrals
    def build_fock_matrix(self, density_matrix):
        fock_matrix = self.h.copy()
        for i in range(self.num_orbitals):
            for j in range(self.num_orbitals):
                fock_matrix[i, j] += np.sum(density_matrix * self.coulomb[i, j])
        return fock_matrix
    def build_density_matrix(self, coefficients):
        density_matrix = np.zeros((self.num_orbitals, self.num_orbitals))
        for i in range(self.num_electrons):
            density_matrix += np.outer(coefficients[:, i], coefficients[:, i])
        return density_matrix
    def diagonalize(self, fock_matrix):
        energy, coefficients = np.linalg.eigh(fock_matrix)
        return energy, coefficients
    def run(self, max_iter=100, tol=1e-6):
        coeffs = np.zeros((self.num_orbitals, self.num_electrons))
        density_matrix = np.zeros((self.num_orbitals, self.num_orbitals))
        for iteration in range(max_iter):
            fock_matrix = self.build_fock_matrix(density_matrix)
            energies, coeffs = self.diagonalize(fock_matrix)
            
            new_density_matrix = self.build_density_matrix(coeffs)
            if np.linalg.norm(new_density_matrix - density_matrix) < tol:
                print(f"Converged in {iteration} iterations")
                break
            density_matrix = new_density_matrix
        return energies, coeffs
# Example usage
if __name__ == "__main__":
    num_electrons = 2
    num_orbitals = 4
    hf = HartreeFock(num_electrons, num_orbitals)
    energies, coefficients = hf.run()
    print("Final energies:", energies)
    print("Final coefficients:", coefficients)

import numpy as np

class SlaterDeterminant:
    def __init__(self, occupied_orbitals=None, n_orbitals=0):
        """
        Initialize a Slater determinant. 
        - `occupied_orbitals` can be an iterable of orbital indices to occupy, or an integer bitstring.
        - `n_orbitals` is the total number of orbitals (useful for bounds and display).
        """
        if occupied_orbitals is None:
            # Default to vacuum state (no occupied orbitals)
            occupied_orbitals = []
        if isinstance(occupied_orbitals, int):
            # If an integer bitstring is given directly
            self.state = occupied_orbitals
        else:
            # If a list/tuple of orbitals is given, set those bits to 1
            self.state = 0
            for orb in occupied_orbitals:
                if orb < 0:
                    raise ValueError("Orbital indices must be non-negative.")
                self.state |= (1 << orb)
        # Determine number of orbitals: if not provided, infer from highest bit
        highest_orbital = self.state.bit_length()  # one more than the index of highest set bit
        self.n_orbitals = max(n_orbitals, highest_orbital)
    
    def create(self, orbital):
        """
        Apply the fermionic creation operator a†_(orbital).
        Returns a tuple (sign, new_determinant) if successful, or None if the result is the zero state.
        """
        if orbital < 0 or orbital >= self.n_orbitals:
            raise IndexError("Orbital index out of range.")
        # If orbital is already occupied, creation yields zero (return None)
        if self.state & (1 << orbital):
            return None
        # Count occupied orbitals with index less than the given orbital (for sign)
        lower_bits_mask = (1 << orbital) - 1       # mask of all bits lower than 'orbital'
        num_occupied_lower = (self.state & lower_bits_mask).bit_count()
        sign = 1 if num_occupied_lower % 2 == 0 else -1
        # Set the orbital bit to 1 to get the new state
        new_state = self.state | (1 << orbital)
        new_det = SlaterDeterminant(new_state, n_orbitals=self.n_orbitals)
        return (sign, new_det)
    
    def annihilate(self, orbital):
        """
        Apply the fermionic annihilation operator a_(orbital).
        Returns a tuple (sign, new_determinant) if successful, or None if the result is the zero state.
        """
        if orbital < 0 or orbital >= self.n_orbitals:
            raise IndexError("Orbital index out of range.")
        # If orbital is not occupied, annihilation yields zero (return None)
        if not (self.state & (1 << orbital)):
            return None
        # Count occupied orbitals with index less than the given orbital (for sign)
        lower_bits_mask = (1 << orbital) - 1       # mask of all bits lower than 'orbital'
        num_occupied_lower = (self.state & lower_bits_mask).bit_count()
        sign = 1 if num_occupied_lower % 2 == 0 else -1
        # Clear the orbital bit to 0 to get the new state
        new_state = self.state & ~(1 << orbital)
        new_det = SlaterDeterminant(new_state, n_orbitals=self.n_orbitals)
        return (sign, new_det)
    
    def __repr__(self):
        # Return a binary string representation padded to n_orbitals bits (orbital 0 is rightmost bit).
        bitstring = format(self.state, '0{}b'.format(self.n_orbitals))
        return f"|{bitstring}>"


# Define a Slater determinant with 5 orbitals, initially occupying orbitals 0, 2, and 4:
det = SlaterDeterminant([0, 2, 4], n_orbitals=5)
print("Initial state:", det)   # should print |10101> (orbitals 4,2,0 occupied)

# Attempt to annihilate a fermion in orbital 1 (which is unoccupied):
result = det.annihilate(1)
print("annihilate(1) result:", result)   # None, since orbital 1 was empty

# Annihilate a fermion in orbital 2 (occupied):
phase, new_det = det.annihilate(2)       # should succeed
print("annihilate(2) -> sign =", phase, ", new state =", new_det)
# Orbital 2 had one occupied orbital below it (orbital 0), so expect sign = -1.
# new_det should be |10001> (orbitals 4 and 0 remain occupied).

# Annihilate a fermion in orbital 4 (occupied):
phase, new_det = det.annihilate(4)
print("annihilate(4) -> sign =", phase, ", new state =", new_det)
# Orbital 4 had two occupied orbitals below (orbitals 0 and 2 were occupied), so sign = +1 (even parity).
# new_det should be |00101> (orbitals 2 and 0 remain after removing orbital 4).

# Create a fermion in orbital 1 (which is currently empty in the original det):
phase, new_det = det.create(1)
print("create(1) -> sign =", phase, ", new state =", new_det)
# Orbital 1 had one occupied orbital below it (orbital 0), so sign = -1.
# new_det should be |11101> (orbitals 4, 2, 1, 0 occupied).

# Create a fermion in orbital 4 (which is already occupied in the original det):
result = det.create(4)
print("create(4) result:", result)      # None, since orbital 4 was already occupied


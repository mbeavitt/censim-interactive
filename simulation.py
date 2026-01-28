"""
Centromere Evolution Simulation Engine

Based on the censim model - simulates evolution of a centromeric repeat array using:
- SNPs (single nucleotide polymorphisms)
- Tandem duplications (duplicate a segment in place)
- Deletions (remove a segment)

All INDEL sizes are frame-aligned (multiples of repeat_size) to maintain
repeat unit boundaries.
"""

import random
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Callable
import copy

# Default CEN178 monomer
DEFAULT_MONOMER = "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG"

BASES = ['A', 'C', 'G', 'T']


class RepeatSequence:
    """Store sequence as a list of repeat units for efficient insertions/deletions.

    Instead of storing millions of characters, we store ~10-15K repeat units.
    Insertions/deletions operate on units, which is ~178x faster.
    """

    def __init__(self, units: List[bytearray], repeat_size: int = 178):
        self.repeat_size = repeat_size
        self.units = units

    @classmethod
    def from_sequence(cls, sequence: str, repeat_size: int = 178) -> 'RepeatSequence':
        """Create from a full DNA sequence string."""
        units = []
        for i in range(0, len(sequence), repeat_size):
            units.append(bytearray(sequence[i:i+repeat_size], 'ascii'))
        return cls(units, repeat_size)

    @classmethod
    def from_monomer(cls, monomer: str, num_copies: int) -> 'RepeatSequence':
        """Create from repeated copies of a monomer."""
        repeat_size = len(monomer)
        units = [bytearray(monomer, 'ascii') for _ in range(num_copies)]
        return cls(units, repeat_size)

    def __len__(self) -> int:
        """Return total sequence length in base pairs."""
        if not self.units:
            return 0
        return len(self.units) * self.repeat_size

    def num_units(self) -> int:
        """Return number of repeat units."""
        return len(self.units)

    def __getitem__(self, idx: int) -> str:
        """Get character at base pair position idx."""
        unit_num = idx // self.repeat_size
        pos_in_unit = idx % self.repeat_size
        return chr(self.units[unit_num][pos_in_unit])

    def __setitem__(self, idx: int, value: str):
        """Set character at base pair position idx."""
        unit_num = idx // self.repeat_size
        pos_in_unit = idx % self.repeat_size
        self.units[unit_num][pos_in_unit] = ord(value)

    def get_unit(self, idx: int) -> str:
        """Get repeat unit at index as string."""
        return self.units[idx].decode('ascii')

    def get_unit_bytes(self, idx: int) -> bytearray:
        """Get repeat unit at index as bytearray."""
        return self.units[idx]

    def duplicate_units(self, start: int, end: int):
        """Tandem duplication: duplicate units[start:end] in place.

        Before: [A][B][C][D][E]
        duplicate_units(1, 3):
        After:  [A][B][C][B][C][D][E]
        """
        # Copy the units to duplicate
        units_to_insert = [bytearray(unit) for unit in self.units[start:end]]
        # Insert after the original segment
        self.units[end:end] = units_to_insert

    def delete_units(self, start: int, end: int):
        """Delete units from start to end.

        Before: [A][B][C][D][E]
        delete_units(1, 3):
        After:  [A][D][E]
        """
        del self.units[start:end]

    def to_string(self) -> str:
        """Convert entire sequence to string."""
        return ''.join(unit.decode('ascii') for unit in self.units)

    def get_unique_sequences(self) -> List[str]:
        """Get list of unique repeat sequences."""
        return list(set(unit.decode('ascii') for unit in self.units))


@dataclass
class SimulationParams:
    """Parameters controlling the simulation."""

    # INDEL rates (expected events per generation, Poisson lambda)
    indel_rate: float = 0.5           # lambda for total INDEL events
    indel_size_lambda: float = 7.6    # lambda for number of repeats per INDEL

    # SNP rate
    snp_rate: float = 0.1             # lambda for SNP events per generation

    # Array bounds
    min_array_size: int = 300         # Minimum number of repeats (collapse threshold)
    max_array_size: int = 50000       # Maximum number of repeats
    bounding_enabled: bool = True     # Whether to enforce bounds

    # Initial state
    initial_size: int = 10000         # Starting number of repeats
    initial_monomer: str = DEFAULT_MONOMER

    # Collapse detection
    max_consecutive_failures: int = 5000  # Max failed attempts before declaring collapse

    def copy(self) -> 'SimulationParams':
        return copy.deepcopy(self)


@dataclass
class MutationEvent:
    """Record of a mutation event."""
    generation: int
    event_type: str  # 'snp', 'dup', 'del'
    position: int    # Unit position (or bp position for SNP)
    size: int = 1    # Number of repeats affected
    details: str = ""


@dataclass
class SimulationState:
    """Current state of the simulation."""
    generation: int = 0
    seq: Optional[RepeatSequence] = None
    history: List[MutationEvent] = field(default_factory=list)
    collapsed: bool = False

    @property
    def size(self) -> int:
        return self.seq.num_units() if self.seq else 0

    @property
    def repeats(self) -> List[str]:
        """Get list of repeat sequences (for colorizer compatibility)."""
        if not self.seq:
            return []
        return [self.seq.get_unit(i) for i in range(self.seq.num_units())]


class CentromereSimulator:
    """
    Simulates centromere repeat array evolution.

    Based on the censim model with frame-aligned INDELs.
    """

    def __init__(self, params: Optional[SimulationParams] = None, seed: int = None):
        self.params = params or SimulationParams()
        self.state = SimulationState()

        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Callbacks for real-time updates
        self.on_mutation: Optional[Callable[[MutationEvent], None]] = None
        self.on_generation: Optional[Callable[[SimulationState], None]] = None

    def initialize(self):
        """Initialize the repeat array with identical monomers."""
        self.state = SimulationState(
            generation=0,
            seq=RepeatSequence.from_monomer(
                self.params.initial_monomer,
                self.params.initial_size
            ),
            history=[],
            collapsed=False
        )

    def _apply_snp(self) -> List[MutationEvent]:
        """Apply SNP mutations for this generation."""
        events = []
        n_snps = np.random.poisson(self.params.snp_rate)

        for _ in range(n_snps):
            if self.state.size == 0:
                break

            # Pick random base pair position
            bp_pos = random.randint(0, len(self.state.seq) - 1)
            old_base = self.state.seq[bp_pos]

            # Choose a different base
            new_base = random.choice([b for b in BASES if b != old_base])
            self.state.seq[bp_pos] = new_base

            event = MutationEvent(
                generation=self.state.generation,
                event_type='snp',
                position=bp_pos,
                size=1,
                details=f"{old_base}->{new_base}"
            )
            events.append(event)

            if self.on_mutation:
                self.on_mutation(event)

        return events

    def _apply_indels(self) -> List[MutationEvent]:
        """Apply INDEL mutations (tandem duplications and deletions)."""
        events = []
        n_indels = np.random.poisson(self.params.indel_rate)
        consecutive_failures = 0

        count = 0
        while count < n_indels:
            if self.state.size == 0:
                self.state.collapsed = True
                break

            # Choose duplication or deletion with equal probability
            is_dup = random.choice([True, False])

            # Sample size in number of repeat units
            indel_size = max(1, int(np.random.poisson(self.params.indel_size_lambda)))

            # Pick random starting position (in repeat units)
            unit_start = random.randint(0, self.state.size - 1)
            unit_end = unit_start + indel_size

            # Bounds check
            if unit_end > self.state.size:
                consecutive_failures += 1
                if consecutive_failures >= self.params.max_consecutive_failures:
                    self.state.collapsed = True
                    break
                continue

            # Check array size limits
            if is_dup:
                # Duplication would grow array
                if self.params.bounding_enabled:
                    if self.state.size + indel_size > self.params.max_array_size:
                        consecutive_failures += 1
                        continue

                self.state.seq.duplicate_units(unit_start, unit_end)
                event = MutationEvent(
                    generation=self.state.generation,
                    event_type='dup',
                    position=unit_start,
                    size=indel_size,
                    details=f"units {unit_start}-{unit_end} duplicated"
                )
            else:
                # Deletion would shrink array
                if self.params.bounding_enabled:
                    if self.state.size - indel_size < self.params.min_array_size:
                        consecutive_failures += 1
                        continue

                self.state.seq.delete_units(unit_start, unit_end)
                event = MutationEvent(
                    generation=self.state.generation,
                    event_type='del',
                    position=unit_start,
                    size=indel_size,
                    details=f"units {unit_start}-{unit_end} deleted"
                )

            # Reset failure counter on success
            consecutive_failures = 0
            events.append(event)
            count += 1

            if self.on_mutation:
                self.on_mutation(event)

        # Check for collapse
        if self.state.size < self.params.min_array_size:
            self.state.collapsed = True

        return events

    def step(self) -> List[MutationEvent]:
        """Advance the simulation by one generation."""
        if self.state.collapsed:
            return []

        self.state.generation += 1
        events = []

        # Apply SNPs
        events.extend(self._apply_snp())

        # Apply INDELs
        events.extend(self._apply_indels())

        # Record history
        self.state.history.extend(events)

        if self.on_generation:
            self.on_generation(self.state)

        return events

    def run(self, generations: int) -> SimulationState:
        """Run simulation for a number of generations."""
        for _ in range(generations):
            if self.state.collapsed:
                break
            self.step()
        return self.state

    def get_unique_sequences(self) -> List[str]:
        """Get list of unique sequences in current array."""
        if not self.state.seq:
            return []
        return self.state.seq.get_unique_sequences()

    def get_statistics(self) -> dict:
        """Get current simulation statistics."""
        unique = self.get_unique_sequences()

        # Count event types
        snps = sum(1 for e in self.state.history if e.event_type == 'snp')
        dups = sum(1 for e in self.state.history if e.event_type == 'dup')
        dels = sum(1 for e in self.state.history if e.event_type == 'del')

        return {
            'generation': self.state.generation,
            'array_size': self.state.size,
            'unique_sequences': len(unique),
            'diversity': len(unique) / self.state.size if self.state.size > 0 else 0,
            'total_events': len(self.state.history),
            'snps': snps,
            'duplications': dups,
            'deletions': dels,
            'collapsed': self.state.collapsed,
        }


# Quick test
if __name__ == "__main__":
    print("Testing CentromereSimulator...")
    print(f"Default monomer: {len(DEFAULT_MONOMER)}bp")

    # Create simulator with smaller initial size for quick test
    params = SimulationParams(initial_size=1000)
    sim = CentromereSimulator(params, seed=42)
    sim.initialize()

    print(f"\nInitial state:")
    print(f"  Array size: {sim.state.size} repeats")
    print(f"  Unique sequences: {len(sim.get_unique_sequences())}")

    # Run for 100 generations
    print(f"\nRunning 100 generations...")
    sim.run(100)

    stats = sim.get_statistics()
    print(f"\nFinal state:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")

    # Show some history
    print(f"\nLast 10 events:")
    for event in sim.state.history[-10:]:
        print(f"  Gen {event.generation}: {event.event_type} at {event.position} (size {event.size}) {event.details}")

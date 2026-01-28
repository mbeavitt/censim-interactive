"""
Centromere Evolution Simulation Engine

Simulates the evolution of a centromeric repeat array over time using:
- SNPs (single nucleotide polymorphisms)
- Insertions (duplication of repeats)
- Deletions (removal of repeats)

All mutation rates are Poisson-distributed and configurable.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Callable
import copy

# Default CEN178 monomer
DEFAULT_MONOMER = "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG"

BASES = ['A', 'C', 'G', 'T']


@dataclass
class SimulationParams:
    """Parameters controlling the simulation."""

    # INDEL rates (per generation)
    insertion_rate: float = 0.25      # lambda for insertion events
    deletion_rate: float = 0.25       # lambda for deletion events
    indel_size_lambda: float = 7.6    # lambda for size of INDELs (Poisson)

    # SNP rate
    snp_rate: float = 0.1             # lambda for SNP events per generation

    # Array bounds
    min_array_size: int = 100         # Minimum number of repeats
    max_array_size: int = 50000       # Maximum number of repeats
    bounding_enabled: bool = True     # Whether to enforce bounds

    # Initial state
    initial_size: int = 10000         # Starting number of repeats
    initial_monomer: str = DEFAULT_MONOMER

    def copy(self) -> 'SimulationParams':
        return copy.deepcopy(self)


@dataclass
class MutationEvent:
    """Record of a mutation event."""
    generation: int
    event_type: str  # 'snp', 'insertion', 'deletion'
    position: int    # Position in array (or base position for SNP)
    size: int = 1    # Number of repeats affected (or 1 for SNP)
    details: str = ""


@dataclass
class SimulationState:
    """Current state of the simulation."""
    generation: int = 0
    repeats: List[str] = field(default_factory=list)
    history: List[MutationEvent] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.repeats)


class CentromereSimulator:
    """
    Simulates centromere repeat array evolution.

    The array is represented as a list of DNA sequences (repeat monomers).
    Each generation, mutations occur according to Poisson-distributed rates.
    """

    def __init__(self, params: Optional[SimulationParams] = None, seed: int = None):
        self.params = params or SimulationParams()
        self.state = SimulationState()
        self.rng = np.random.default_rng(seed)

        # Callbacks for real-time updates
        self.on_mutation: Optional[Callable[[MutationEvent], None]] = None
        self.on_generation: Optional[Callable[[SimulationState], None]] = None

    def initialize(self):
        """Initialize the repeat array with identical monomers."""
        self.state = SimulationState(
            generation=0,
            repeats=[self.params.initial_monomer] * self.params.initial_size,
            history=[]
        )

    def _apply_snp(self, repeat_idx: int) -> MutationEvent:
        """Apply a single nucleotide polymorphism to a repeat."""
        seq = list(self.state.repeats[repeat_idx])
        pos = self.rng.integers(0, len(seq))
        old_base = seq[pos]

        # Choose a different base
        new_base = self.rng.choice([b for b in BASES if b != old_base])
        seq[pos] = new_base

        self.state.repeats[repeat_idx] = ''.join(seq)

        return MutationEvent(
            generation=self.state.generation,
            event_type='snp',
            position=repeat_idx,
            size=1,
            details=f"pos {pos}: {old_base}->{new_base}"
        )

    def _apply_insertion(self, position: int, size: int) -> Optional[MutationEvent]:
        """Insert duplicated repeats at a position."""
        if self.params.bounding_enabled:
            if self.state.size + size > self.params.max_array_size:
                size = self.params.max_array_size - self.state.size
                if size <= 0:
                    return None

        # Duplicate repeats from a nearby region
        source_start = max(0, position - size)
        source_end = min(self.state.size, source_start + size)
        to_insert = self.state.repeats[source_start:source_end]

        # If we couldn't get enough, pad with copies
        while len(to_insert) < size:
            to_insert.append(self.state.repeats[position % self.state.size])

        # Insert at position
        self.state.repeats = (
            self.state.repeats[:position] +
            to_insert[:size] +
            self.state.repeats[position:]
        )

        return MutationEvent(
            generation=self.state.generation,
            event_type='insertion',
            position=position,
            size=size,
            details=f"duplicated from {source_start}"
        )

    def _apply_deletion(self, position: int, size: int) -> Optional[MutationEvent]:
        """Delete repeats starting at a position."""
        if self.params.bounding_enabled:
            if self.state.size - size < self.params.min_array_size:
                size = self.state.size - self.params.min_array_size
                if size <= 0:
                    return None

        # Ensure we don't delete past the end
        size = min(size, self.state.size - position)
        if size <= 0:
            return None

        self.state.repeats = (
            self.state.repeats[:position] +
            self.state.repeats[position + size:]
        )

        return MutationEvent(
            generation=self.state.generation,
            event_type='deletion',
            position=position,
            size=size,
            details=""
        )

    def step(self) -> List[MutationEvent]:
        """
        Advance the simulation by one generation.

        Returns list of mutation events that occurred.
        """
        events = []

        # Sample number of each event type from Poisson
        n_insertions = self.rng.poisson(self.params.insertion_rate)
        n_deletions = self.rng.poisson(self.params.deletion_rate)
        n_snps = self.rng.poisson(self.params.snp_rate)

        # Apply insertions
        for _ in range(n_insertions):
            if self.state.size == 0:
                break
            pos = self.rng.integers(0, self.state.size)
            size = max(1, self.rng.poisson(self.params.indel_size_lambda))
            event = self._apply_insertion(pos, size)
            if event:
                events.append(event)
                if self.on_mutation:
                    self.on_mutation(event)

        # Apply deletions
        for _ in range(n_deletions):
            if self.state.size <= self.params.min_array_size:
                break
            pos = self.rng.integers(0, self.state.size)
            size = max(1, self.rng.poisson(self.params.indel_size_lambda))
            event = self._apply_deletion(pos, size)
            if event:
                events.append(event)
                if self.on_mutation:
                    self.on_mutation(event)

        # Apply SNPs
        for _ in range(n_snps):
            if self.state.size == 0:
                break
            repeat_idx = self.rng.integers(0, self.state.size)
            event = self._apply_snp(repeat_idx)
            events.append(event)
            if self.on_mutation:
                self.on_mutation(event)

        # Update state
        self.state.generation += 1
        self.state.history.extend(events)

        if self.on_generation:
            self.on_generation(self.state)

        return events

    def run(self, generations: int) -> SimulationState:
        """Run simulation for a number of generations."""
        for _ in range(generations):
            self.step()
        return self.state

    def get_unique_sequences(self) -> List[str]:
        """Get list of unique sequences in current array."""
        return list(set(self.state.repeats))

    def get_statistics(self) -> dict:
        """Get current simulation statistics."""
        unique = self.get_unique_sequences()
        return {
            'generation': self.state.generation,
            'array_size': self.state.size,
            'unique_sequences': len(unique),
            'diversity': len(unique) / self.state.size if self.state.size > 0 else 0,
            'total_events': len(self.state.history),
            'snps': sum(1 for e in self.state.history if e.event_type == 'snp'),
            'insertions': sum(1 for e in self.state.history if e.event_type == 'insertion'),
            'deletions': sum(1 for e in self.state.history if e.event_type == 'deletion'),
        }


# Quick test
if __name__ == "__main__":
    print("Testing CentromereSimulator...")

    # Create simulator with smaller initial size for quick test
    params = SimulationParams(initial_size=100)
    sim = CentromereSimulator(params, seed=42)
    sim.initialize()

    print(f"\nInitial state:")
    print(f"  Array size: {sim.state.size}")
    print(f"  Unique sequences: {len(sim.get_unique_sequences())}")

    # Run for 100 generations
    print(f"\nRunning 100 generations...")
    sim.run(100)

    stats = sim.get_statistics()
    print(f"\nFinal state:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    # Show some history
    print(f"\nLast 5 events:")
    for event in sim.state.history[-5:]:
        print(f"  Gen {event.generation}: {event.event_type} at pos {event.position} (size {event.size})")

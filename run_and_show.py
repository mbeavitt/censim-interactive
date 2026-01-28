#!/usr/bin/env python3
"""
Run a long simulation and show the final state.
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt

from simulation import CentromereSimulator, SimulationParams, DEFAULT_MONOMER
from colorizer import OrthogonalProjectionColorizer


def run_simulation(generations: int, seed: int = 42) -> CentromereSimulator:
    """Run simulation for specified generations."""
    print(f"Running simulation for {generations:,} generations...")

    params = SimulationParams(initial_size=10000)
    sim = CentromereSimulator(params, seed=seed)
    sim.initialize()

    start = time.time()

    # Run in chunks and report progress
    chunk_size = min(100_000, generations // 10) if generations > 100_000 else generations
    completed = 0

    while completed < generations:
        to_run = min(chunk_size, generations - completed)
        sim.run(to_run)
        completed += to_run

        stats = sim.get_statistics()
        elapsed = time.time() - start

        if generations > 100_000:
            print(f"  Gen {stats['generation']:,}: size={stats['array_size']:,}, "
                  f"unique={stats['unique_sequences']:,}, "
                  f"diversity={stats['diversity']:.4f} ({elapsed:.1f}s)")

        if stats['collapsed']:
            print("  Array collapsed!")
            break

    print(f"\nSimulation complete in {time.time() - start:.1f}s")
    return sim


def visualize_result(sim: CentromereSimulator, grid_width: int = 100):
    """Visualize the final simulation state."""
    print("\nGenerating visualization...")

    # Create colorizer
    colorizer = OrthogonalProjectionColorizer(seq_len=len(DEFAULT_MONOMER))

    # Get unique sequences and compute colors
    unique_seqs = sim.get_unique_sequences()
    colorizer.colorize_batch(unique_seqs)

    # Build RGB array
    n_repeats = sim.state.size
    n_rows = (n_repeats + grid_width - 1) // grid_width

    rgb = np.ones((n_rows, grid_width, 3)) * 0.15  # Dark background

    for idx, seq in enumerate(sim.state.repeats):
        row = idx // grid_width
        col = idx % grid_width
        rgb[row, col] = colorizer.get_color(seq)

    # Create figure
    stats = sim.get_statistics()

    fig, (ax_main, ax_stats) = plt.subplots(1, 2, figsize=(18, 10),
                                             gridspec_kw={'width_ratios': [3, 1]})

    # Main visualization
    ax_main.imshow(rgb, aspect='auto', interpolation='nearest')
    ax_main.set_title(f"Centromere Evolution - {stats['generation']:,} Generations\n"
                      f"{stats['array_size']:,} repeats, {stats['unique_sequences']:,} unique",
                      fontsize=14)
    ax_main.axis('off')

    # Stats panel
    ax_stats.axis('off')
    stats_text = (
        f"Simulation Results\n"
        f"{'=' * 30}\n\n"
        f"Generations: {stats['generation']:,}\n"
        f"Array size: {stats['array_size']:,}\n"
        f"Unique sequences: {stats['unique_sequences']:,}\n"
        f"Diversity: {stats['diversity']:.4f}\n\n"
        f"Mutation Events\n"
        f"{'-' * 30}\n"
        f"SNPs: {stats['snps']:,}\n"
        f"Duplications: {stats['duplications']:,}\n"
        f"Deletions: {stats['deletions']:,}\n"
        f"Total: {stats['total_events']:,}\n\n"
        f"Array collapsed: {stats['collapsed']}"
    )
    ax_stats.text(0.1, 0.9, stats_text, transform=ax_stats.transAxes,
                  fontsize=12, verticalalignment='top', fontfamily='monospace',
                  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig('simulation_result.png', dpi=150, bbox_inches='tight')
    print(f"Saved to simulation_result.png")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Run simulation and visualize result')
    parser.add_argument('-g', '--generations', type=int, default=2_000_000,
                        help='Number of generations to simulate (default: 2,000,000)')
    parser.add_argument('-s', '--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('-w', '--grid-width', type=int, default=100,
                        help='Grid width for visualization (default: 100)')

    args = parser.parse_args()

    sim = run_simulation(args.generations, args.seed)
    visualize_result(sim, args.grid_width)


if __name__ == "__main__":
    main()

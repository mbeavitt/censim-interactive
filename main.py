#!/usr/bin/env python3
"""
Centromere Evolution Simulator - Main Entry Point

A real-time visualization of centromeric repeat array evolution.
This is the Python prototype before the C/raylib implementation.

Usage:
    python main.py              # Run interactive visualizer
    python main.py --headless   # Run simulation without GUI (for testing)
    python main.py --benchmark  # Benchmark simulation speed
"""

import argparse
import time


def run_visualizer():
    """Run the interactive matplotlib visualizer."""
    from visualizer import SimulationVisualizer

    print("=" * 60)
    print("Centromere Evolution Simulator")
    print("=" * 60)
    print()
    print("Python prototype - will be reimplemented in C with raylib")
    print()

    viz = SimulationVisualizer(grid_width=100)
    viz.show()


def run_headless(generations: int = 1000):
    """Run simulation without visualization."""
    from simulation import CentromereSimulator, SimulationParams

    print("Running headless simulation...")

    params = SimulationParams(initial_size=10000)
    sim = CentromereSimulator(params, seed=42)
    sim.initialize()

    print(f"Initial: {sim.state.size:,} repeats")

    start = time.time()
    for i in range(generations):
        sim.step()
        if (i + 1) % 100 == 0:
            stats = sim.get_statistics()
            print(f"Gen {i+1:4d}: size={stats['array_size']:,}, "
                  f"unique={stats['unique_sequences']:,}, "
                  f"diversity={stats['diversity']:.3f}")

    elapsed = time.time() - start
    stats = sim.get_statistics()

    print()
    print(f"Final state after {generations} generations:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value:,}" if isinstance(value, int) else f"  {key}: {value}")

    print(f"\nElapsed time: {elapsed:.2f}s ({generations/elapsed:.1f} gen/s)")


def run_benchmark():
    """Benchmark simulation and colorization speed."""
    from simulation import CentromereSimulator, SimulationParams
    from colorizer import OrthogonalProjectionColorizer, DEFAULT_MONOMER

    print("=" * 60)
    print("Benchmark")
    print("=" * 60)

    # Simulation benchmark
    print("\n1. Simulation speed (no colorization)")
    params = SimulationParams(initial_size=10000)
    sim = CentromereSimulator(params, seed=42)
    sim.initialize()

    n_gens = 1000
    start = time.time()
    sim.run(n_gens)
    elapsed = time.time() - start

    print(f"   {n_gens} generations in {elapsed:.2f}s")
    print(f"   {n_gens/elapsed:.1f} generations/second")
    print(f"   Final size: {sim.state.size:,}")

    # Colorization benchmark
    print("\n2. Colorization speed")
    unique_seqs = sim.get_unique_sequences()
    colorizer = OrthogonalProjectionColorizer(seq_len=len(DEFAULT_MONOMER))

    start = time.time()
    colorizer.colorize_batch(unique_seqs)
    elapsed = time.time() - start

    print(f"   {len(unique_seqs):,} unique sequences in {elapsed:.3f}s")
    print(f"   {len(unique_seqs)/elapsed:.0f} sequences/second")

    # Per-sequence lookup benchmark
    print("\n3. Individual color lookup speed")
    n_lookups = 100000
    start = time.time()
    for i in range(n_lookups):
        seq = sim.state.repeats[i % len(sim.state.repeats)]
        _ = colorizer.get_color(seq)
    elapsed = time.time() - start

    print(f"   {n_lookups:,} lookups in {elapsed:.3f}s")
    print(f"   {n_lookups/elapsed:.0f} lookups/second")

    # Combined benchmark (simulate what real-time viz needs)
    print("\n4. Combined (simulate + colorize full array)")
    sim2 = CentromereSimulator(SimulationParams(initial_size=10000), seed=123)
    sim2.initialize()
    colorizer2 = OrthogonalProjectionColorizer(seq_len=len(DEFAULT_MONOMER))

    n_frames = 100
    start = time.time()
    for _ in range(n_frames):
        sim2.step()
        # Color all repeats (what visualization needs)
        for seq in sim2.state.repeats:
            _ = colorizer2.get_color(seq)
    elapsed = time.time() - start

    print(f"   {n_frames} frames in {elapsed:.2f}s")
    print(f"   {n_frames/elapsed:.1f} FPS (frames/second)")
    print(f"   Target for raylib: 60 FPS")


def main():
    parser = argparse.ArgumentParser(
        description='Centromere Evolution Simulator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This is the Python prototype for testing the simulation logic.
The final version will be implemented in C with raylib for real-time performance.
        """
    )

    parser.add_argument('--headless', action='store_true',
                        help='Run simulation without GUI')
    parser.add_argument('--benchmark', action='store_true',
                        help='Run performance benchmarks')
    parser.add_argument('--generations', '-g', type=int, default=1000,
                        help='Number of generations for headless mode')

    args = parser.parse_args()

    if args.benchmark:
        run_benchmark()
    elif args.headless:
        run_headless(args.generations)
    else:
        run_visualizer()


if __name__ == "__main__":
    main()

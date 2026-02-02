#!/usr/bin/env python3
"""
UMAP-based visualization for centromere simulation output.

Takes a FASTA file of repeat sequences and creates a 2D grid visualization
with colors derived from UMAP projection of Levenshtein distance matrix.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.preprocessing import MinMaxScaler
from umap import UMAP
import Levenshtein
import warnings


def read_fasta(fasta_path):
    """Read sequences from FASTA file."""
    sequences = []
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append(''.join(current_seq))

    return sequences


def compute_color_map(unique_sequences):
    """Compute UMAP-based color mapping for unique sequences."""
    n_unique = len(unique_sequences)
    print(f"  Computing distance matrix ({n_unique} x {n_unique})...")

    # Build distance matrix
    distance_matrix = np.zeros((n_unique, n_unique))
    for i in range(n_unique):
        for j in range(i + 1, n_unique):
            dist = Levenshtein.distance(unique_sequences[i], unique_sequences[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

        if (i + 1) % 100 == 0:
            print(f"    {i + 1}/{n_unique} rows...")

    # Project to 3D RGB space using UMAP
    print(f"  Projecting to 3D color space using UMAP...")
    warnings.filterwarnings('ignore')

    reducer = UMAP(
        n_components=3,
        metric='precomputed',
        n_neighbors=min(15, n_unique - 1),
        min_dist=0.1,
        n_jobs=1,
        init='random',
        random_state=42
    )
    rgb_coords = reducer.fit_transform(distance_matrix)

    # Normalize to [0, 1]
    scaler = MinMaxScaler()
    rgb_coords = scaler.fit_transform(rgb_coords)
    rgb_coords = np.clip(rgb_coords, 0, 1)

    return {seq: rgb_coords[i] for i, seq in enumerate(unique_sequences)}, rgb_coords


def save_3d_plot(rgb_coords, unique_sequences, output_dir):
    """Save 3D visualization of color space."""
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(
        rgb_coords[:, 0], rgb_coords[:, 1], rgb_coords[:, 2],
        c=rgb_coords, s=50, alpha=0.8, edgecolors='k', linewidth=0.5
    )

    ax.set_xlabel('Red (UMAP 1)')
    ax.set_ylabel('Green (UMAP 2)')
    ax.set_zlabel('Blue (UMAP 3)')
    ax.set_title(f'UMAP Color Space: {len(unique_sequences)} Unique Sequences')

    output_path = output_dir / 'umap_3d_color_space.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def visualize_2d(sequences, color_map, output_path=None, grid_width=100, tile_size=8, show=False):
    """Create 2D grid visualization."""
    n_repeats = len(sequences)
    n_rows = (n_repeats + grid_width - 1) // grid_width

    print(f"  Creating {grid_width}x{n_rows} grid...")

    # Build RGB array
    height_px = n_rows * tile_size
    width_px = grid_width * tile_size
    rgb_array = np.ones((height_px, width_px, 3)) * 0.15  # Dark background

    for idx, seq in enumerate(sequences):
        row = idx // grid_width
        col = idx % grid_width

        y_start = row * tile_size
        y_end = (row + 1) * tile_size
        x_start = col * tile_size
        x_end = (col + 1) * tile_size

        if seq in color_map:
            rgb_array[y_start:y_end, x_start:x_end] = color_map[seq]

    # Create figure
    fig_width = max(12, grid_width * tile_size / 100)
    fig_height = max(8, n_rows * tile_size / 100)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.imshow(rgb_array, aspect='auto', interpolation='nearest')
    ax.axis('off')

    title = f'Simulation Visualization (UMAP colors)\n{n_repeats:,} repeats in {grid_width}x{n_rows} grid'
    plt.suptitle(title, fontsize=14, y=0.98)
    plt.tight_layout()

    if show:
        print("  Showing plot...")
        # Signal that we're ready to show (for UI feedback)
        signal_file = Path("/tmp/censim_umap_ready")
        signal_file.touch()
        plt.show()
    elif output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='UMAP visualization for simulation output')
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output image path (if not using --show)')
    parser.add_argument('-w', '--grid-width', type=int, default=100, help='Grid width')
    parser.add_argument('-t', '--tile-size', type=int, default=8, help='Tile size in pixels')
    parser.add_argument('--no-3d', action='store_true', help='Skip 3D color space plot')
    parser.add_argument('--show', action='store_true', help='Show plot interactively instead of saving')

    args = parser.parse_args()

    print(f"Loading {args.fasta}...")
    sequences = read_fasta(args.fasta)
    print(f"  Loaded {len(sequences):,} sequences")

    unique_sequences = sorted(set(sequences))
    print(f"  {len(unique_sequences):,} unique sequences")

    print("Computing UMAP color mapping...")
    color_map, rgb_coords = compute_color_map(unique_sequences)

    output_path = Path(args.output) if args.output else None
    if output_path:
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)

        if not args.no_3d:
            save_3d_plot(rgb_coords, unique_sequences, output_dir)

    print("Creating 2D visualization...")
    visualize_2d(sequences, color_map, output_path, args.grid_width, args.tile_size, show=args.show)

    print("Done!")


if __name__ == '__main__':
    main()

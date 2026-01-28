# Centromere Evolution Simulator

Real-time visualization of centromeric repeat array evolution.

## Overview

This project simulates how centromeric DNA repeat arrays evolve over time through:
- **SNPs** (single nucleotide polymorphisms)
- **Insertions** (tandem duplications)
- **Deletions** (repeat removal)

All mutation rates follow Poisson distributions with configurable lambda parameters.

## Project Structure

```
censim-interactive/
├── Python prototype (current)
│   ├── main.py          # Entry point
│   ├── simulation.py    # Evolution engine
│   ├── colorizer.py     # DNA sequence to RGB color
│   └── visualizer.py    # Matplotlib-based visualization
│
└── C/raylib version (planned)
    └── ...
```

## Python Prototype

### Requirements

```bash
pip install numpy matplotlib
```

### Usage

```bash
# Interactive visualizer
python main.py

# Headless simulation (for testing)
python main.py --headless -g 1000

# Performance benchmark
python main.py --benchmark
```

### Controls (Interactive Mode)

- **Start/Stop**: Toggle simulation running
- **Step**: Advance one generation
- **Reset**: Restart with initial conditions
- **Sliders**: Adjust mutation rates in real-time
  - Insertion rate (lambda)
  - Deletion rate (lambda)
  - INDEL size (lambda for Poisson-distributed size)
  - SNP rate (lambda)
  - Animation speed
- **Bounding**: Toggle min/max array size enforcement

## Simulation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `insertion_rate` | 0.25 | Expected insertions per generation |
| `deletion_rate` | 0.25 | Expected deletions per generation |
| `indel_size_lambda` | 7.6 | Expected size of each INDEL |
| `snp_rate` | 0.1 | Expected SNPs per generation |
| `min_array_size` | 100 | Minimum repeats (when bounding enabled) |
| `max_array_size` | 50000 | Maximum repeats (when bounding enabled) |
| `initial_size` | 10000 | Starting number of repeats |

## Coloring Method

Uses **orthogonal random projection** to map DNA sequences to RGB colors:

1. One-hot encode each sequence (178bp × 4 bases = 712 dimensions)
2. Project through orthogonalized random matrix (712 → 3 dimensions)
3. Normalize to [0, 1] RGB range

Similar sequences get similar colors, making patterns visible in the visualization.

## Default Monomer

The simulation starts with copies of the CEN178 consensus sequence:
```
AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAA
GACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGA
AGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG
```

## Future: C/raylib Version

The Python prototype will be reimplemented in C with raylib for:
- 60+ FPS real-time rendering
- Larger arrays (100K+ repeats)
- GPU-accelerated colorization
- Interactive parameter adjustment

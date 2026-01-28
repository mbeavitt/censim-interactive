# Centromere Evolution Simulator

Real-time visualization of centromeric repeat array evolution using raylib.

## Overview

Simulates how centromeric DNA repeat arrays evolve over time through:
- **SNPs** (single nucleotide polymorphisms)
- **Tandem duplications** (duplicate a segment in place)
- **Deletions** (remove a segment)

All mutation rates follow Poisson distributions. INDELs are frame-aligned to maintain repeat unit boundaries.

## Building

Requires raylib 5.x and premake5:

```bash
brew install raylib premake

# Generate makefiles and build
premake5 gmake2
make -C build config=release

# Run
./bin/Release/censim
```

## Controls

- **F11 / F**: Toggle fullscreen
- **Start/Stop**: Toggle simulation
- **Step 1000**: Advance 1000 generations
- **Reset**: Restart simulation
- **Sliders**: Adjust mutation rates in real-time
- **Bounding**: Toggle min/max array size enforcement

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `indel_rate` | 0.5 | Expected INDELs per generation |
| `indel_size_lambda` | 7.6 | Expected size of each INDEL |
| `snp_rate` | 0.1 | Expected SNPs per generation |
| `min_array_size` | 300 | Minimum repeats (when bounding enabled) |
| `max_array_size` | 50000 | Maximum repeats |
| `initial_size` | 10000 | Starting number of repeats |

## Coloring Method

Uses **orthogonal random projection** to map DNA sequences to RGB colors:

1. One-hot encode each sequence (178bp × 4 bases = 712 dimensions)
2. Project through orthogonalized random matrix (712 → 3 dimensions)
3. Normalize to [0, 1] RGB range

Similar sequences get similar colors, making evolutionary patterns visible.

## Default Monomer

The simulation starts with copies of the CEN178 consensus sequence (178bp).

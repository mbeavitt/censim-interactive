"""
Orthogonal Random Projection Colorizer for DNA sequences.

Maps DNA sequences to RGB colors using one-hot encoding projected through
an orthogonalized random matrix. Similar sequences get similar colors.
"""

import numpy as np
from typing import Dict, List

# Default CEN178 monomer sequence
DEFAULT_MONOMER = "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG"


class OrthogonalProjectionColorizer:
    """
    Color sequences using orthogonalized random projection from one-hot encoding.

    Uses QR decomposition to make the 3 projection directions orthogonal,
    giving good color spread across RGB space.
    """

    def __init__(self, seq_len: int = 178, seed: int = 42, fixed_bounds: bool = True):
        self.seq_len = seq_len
        self.seed = seed
        self.fixed_bounds = fixed_bounds

        # Create and cache projection matrix
        np.random.seed(seed)
        random_matrix = np.random.randn(seq_len * 4, 3)
        q, _ = np.linalg.qr(random_matrix)
        self.projection = q  # Orthonormal columns

        self.base_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        # Cache for normalized colors
        self._color_cache: Dict[str, np.ndarray] = {}
        self._raw_colors: Dict[str, np.ndarray] = {}

        # Fixed normalization bounds (computed from projection properties)
        # For orthonormal projection, output range is roughly [-sqrt(n), sqrt(n)]
        # where n = number of 1s in one-hot encoding = seq_len
        if fixed_bounds:
            # Estimate bounds empirically by sampling random sequences
            self._compute_fixed_bounds()
        else:
            self._min_vals = None
            self._max_vals = None

    def _compute_fixed_bounds(self):
        """Compute fixed normalization bounds by sampling random sequences."""
        np.random.seed(self.seed + 1000)  # Different seed for sampling

        # Generate many random sequences to estimate the range
        n_samples = 1000
        raw_values = []

        for _ in range(n_samples):
            # Random sequence
            seq = ''.join(np.random.choice(list('ACGT'), self.seq_len))
            one_hot = self._one_hot(seq)
            rgb = one_hot @ self.projection
            raw_values.append(rgb)

        raw_values = np.array(raw_values)

        # Use percentiles to avoid outliers, with some padding
        self._min_vals = np.percentile(raw_values, 1, axis=0) - 0.5
        self._max_vals = np.percentile(raw_values, 99, axis=0) + 0.5

    def _one_hot(self, seq: str) -> np.ndarray:
        """Convert sequence to one-hot encoding."""
        enc = np.zeros(len(seq) * 4)
        for i, base in enumerate(seq):
            if base in self.base_idx:
                enc[i * 4 + self.base_idx[base]] = 1
        return enc

    def _get_projection(self, seq_len: int) -> np.ndarray:
        """Get projection matrix for given sequence length."""
        if seq_len == self.seq_len:
            return self.projection

        # Generate new orthogonal projection for different length
        np.random.seed(self.seed)
        random_matrix = np.random.randn(seq_len * 4, 3)
        q, _ = np.linalg.qr(random_matrix)
        return q

    def get_raw_color(self, seq: str) -> np.ndarray:
        """Get raw (unnormalized) RGB values for a sequence."""
        if seq in self._raw_colors:
            return self._raw_colors[seq]

        one_hot = self._one_hot(seq)
        projection = self._get_projection(len(seq))
        rgb = one_hot @ projection

        self._raw_colors[seq] = rgb
        return rgb

    def colorize_batch(self, sequences: List[str]) -> Dict[str, np.ndarray]:
        """
        Compute normalized RGB colors for a batch of sequences.

        With fixed_bounds=True (default), uses pre-computed bounds for consistent colors.
        With fixed_bounds=False, recomputes bounds from this batch.
        """
        # Compute raw colors
        raw_colors = np.array([self.get_raw_color(seq) for seq in sequences])

        # Update normalization bounds only if not using fixed bounds
        if not self.fixed_bounds:
            self._min_vals = raw_colors.min(axis=0)
            self._max_vals = raw_colors.max(axis=0)

        # Normalize using current bounds
        range_vals = self._max_vals - self._min_vals
        range_vals[range_vals == 0] = 1  # Avoid division by zero

        normalized = (raw_colors - self._min_vals) / range_vals
        normalized = np.clip(normalized, 0, 1)

        # Cache results
        for i, seq in enumerate(sequences):
            self._color_cache[seq] = normalized[i]

        return {seq: normalized[i] for i, seq in enumerate(sequences)}

    def get_color(self, seq: str) -> np.ndarray:
        """
        Get normalized RGB color for a single sequence.

        If the sequence was seen in colorize_batch(), returns cached result.
        Otherwise, computes and caches the color using fixed bounds.
        """
        if seq in self._color_cache:
            return self._color_cache[seq]

        raw = self.get_raw_color(seq)

        # Normalize using fixed bounds
        range_vals = self._max_vals - self._min_vals
        range_vals[range_vals == 0] = 1

        normalized = (raw - self._min_vals) / range_vals
        normalized = np.clip(normalized, 0, 1)

        self._color_cache[seq] = normalized
        return normalized

    def get_color_uint8(self, seq: str) -> tuple:
        """Get RGB color as 0-255 integers (for graphics APIs)."""
        rgb = self.get_color(seq)
        return tuple((rgb * 255).astype(np.uint8))

    def clear_cache(self):
        """Clear the color cache (keeps fixed bounds)."""
        self._color_cache.clear()
        self._raw_colors.clear()
        # Note: does NOT clear _min_vals/_max_vals so colors stay consistent


# Quick test
if __name__ == "__main__":
    colorizer = OrthogonalProjectionColorizer(seq_len=len(DEFAULT_MONOMER))

    # Test with default monomer and some mutations
    test_seqs = [DEFAULT_MONOMER]

    # Add some mutated versions
    for i in range(5):
        mutated = list(DEFAULT_MONOMER)
        mutated[i * 10] = 'T' if mutated[i * 10] != 'T' else 'A'
        test_seqs.append(''.join(mutated))

    colors = colorizer.colorize_batch(test_seqs)

    print(f"Default monomer length: {len(DEFAULT_MONOMER)}")
    print(f"\nColors for test sequences:")
    for i, seq in enumerate(test_seqs):
        rgb = colors[seq]
        rgb_int = colorizer.get_color_uint8(seq)
        print(f"  Seq {i}: RGB = [{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}] = {rgb_int}")

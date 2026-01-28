"""
Prototype Visualizer using Matplotlib

Displays the repeat array as a 2D grid of colored tiles, similar to
the centromere 2D visualization. This is a prototype before moving to raylib.

Uses matplotlib's animation for real-time updates (slow but functional).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.widgets import Slider, Button, CheckButtons
from typing import Optional
import time

from colorizer import OrthogonalProjectionColorizer, DEFAULT_MONOMER
from simulation import CentromereSimulator, SimulationParams


class SimulationVisualizer:
    """
    Matplotlib-based visualizer for the centromere simulation.

    Displays repeats as a 2D grid of colored tiles with interactive controls.
    """

    def __init__(self, grid_width: int = 100, tile_size: int = 4):
        self.grid_width = grid_width
        self.tile_size = tile_size

        # Simulation components
        self.params = SimulationParams(initial_size=10000)
        self.sim: Optional[CentromereSimulator] = None
        self.colorizer = OrthogonalProjectionColorizer(seq_len=len(DEFAULT_MONOMER))

        # State
        self.running = False
        self.generation_delay = 0.05  # seconds between generations

        # Setup will be called when show() is invoked
        self.fig = None
        self.ax_main = None
        self.ax_stats = None
        self.image = None

    def _build_rgb_array(self) -> np.ndarray:
        """Build RGB array from current simulation state."""
        n_repeats = self.sim.state.size
        n_rows = (n_repeats + self.grid_width - 1) // self.grid_width

        # Create RGB array
        rgb = np.ones((n_rows, self.grid_width, 3)) * 0.2  # Dark gray background

        for idx, seq in enumerate(self.sim.state.repeats):
            row = idx // self.grid_width
            col = idx % self.grid_width
            rgb[row, col] = self.colorizer.get_color(seq)

        return rgb

    def _update_display(self):
        """Update the visualization."""
        if self.sim is None:
            return

        # Update main image
        rgb = self._build_rgb_array()
        self.image.set_data(rgb)
        self.image.set_extent([0, self.grid_width, rgb.shape[0], 0])
        self.ax_main.set_ylim(rgb.shape[0], 0)

        # Update stats text
        stats = self.sim.get_statistics()
        stats_text = (
            f"Generation: {stats['generation']}\n"
            f"Array size: {stats['array_size']:,}\n"
            f"Unique seqs: {stats['unique_sequences']:,}\n"
            f"Diversity: {stats['diversity']:.3f}\n"
            f"SNPs: {stats['snps']}\n"
            f"Duplications: {stats['duplications']}\n"
            f"Deletions: {stats['deletions']}\n"
            f"Collapsed: {stats['collapsed']}"
        )
        self.stats_text.set_text(stats_text)

        # Update title
        self.ax_main.set_title(
            f"Centromere Evolution - Gen {stats['generation']} - "
            f"{stats['array_size']:,} repeats"
        )

    def _setup_figure(self):
        """Set up the matplotlib figure with controls."""
        self.fig = plt.figure(figsize=(16, 10))

        # Main visualization area
        self.ax_main = self.fig.add_axes([0.05, 0.25, 0.65, 0.70])
        self.ax_main.set_title("Centromere Evolution Simulation")
        self.ax_main.axis('off')

        # Stats panel
        self.ax_stats = self.fig.add_axes([0.75, 0.55, 0.22, 0.40])
        self.ax_stats.axis('off')
        self.stats_text = self.ax_stats.text(
            0.05, 0.95, "", transform=self.ax_stats.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        )

        # Initialize empty image
        rgb = np.ones((100, self.grid_width, 3)) * 0.2
        self.image = self.ax_main.imshow(rgb, aspect='auto', interpolation='nearest')

        # Sliders
        slider_color = 'lightgoldenrodyellow'

        # INDEL rate slider (combined dup + del)
        ax_indel = self.fig.add_axes([0.75, 0.45, 0.20, 0.03], facecolor=slider_color)
        self.slider_indel = Slider(ax_indel, 'INDEL rate', 0, 3, valinit=self.params.indel_rate)
        self.slider_indel.on_changed(self._on_indel_change)

        # INDEL size slider
        ax_size = self.fig.add_axes([0.75, 0.40, 0.20, 0.03], facecolor=slider_color)
        self.slider_size = Slider(ax_size, 'INDEL size', 1, 30, valinit=self.params.indel_size_lambda)
        self.slider_size.on_changed(self._on_size_change)

        # SNP rate slider
        ax_snp = self.fig.add_axes([0.75, 0.35, 0.20, 0.03], facecolor=slider_color)
        self.slider_snp = Slider(ax_snp, 'SNP rate', 0, 1, valinit=self.params.snp_rate)
        self.slider_snp.on_changed(self._on_snp_change)

        # Speed slider
        ax_speed = self.fig.add_axes([0.75, 0.30, 0.20, 0.03], facecolor=slider_color)
        self.slider_speed = Slider(ax_speed, 'Speed', 0.01, 0.5, valinit=self.generation_delay)
        self.slider_speed.on_changed(self._on_speed_change)

        # Buttons
        ax_start = self.fig.add_axes([0.75, 0.20, 0.10, 0.05])
        self.btn_start = Button(ax_start, 'Start/Stop')
        self.btn_start.on_clicked(self._on_start_stop)

        ax_reset = self.fig.add_axes([0.86, 0.20, 0.10, 0.05])
        self.btn_reset = Button(ax_reset, 'Reset')
        self.btn_reset.on_clicked(self._on_reset)

        ax_step = self.fig.add_axes([0.75, 0.13, 0.10, 0.05])
        self.btn_step = Button(ax_step, 'Step')
        self.btn_step.on_clicked(self._on_step)

        # Checkbox for bounding
        ax_check = self.fig.add_axes([0.86, 0.10, 0.12, 0.10])
        self.check_bound = CheckButtons(ax_check, ['Bounding'], [self.params.bounding_enabled])
        self.check_bound.on_clicked(self._on_bound_toggle)

    def _on_indel_change(self, val):
        self.params.indel_rate = val
        if self.sim:
            self.sim.params.indel_rate = val

    def _on_size_change(self, val):
        self.params.indel_size_lambda = val
        if self.sim:
            self.sim.params.indel_size_lambda = val

    def _on_snp_change(self, val):
        self.params.snp_rate = val
        if self.sim:
            self.sim.params.snp_rate = val

    def _on_speed_change(self, val):
        self.generation_delay = val

    def _on_start_stop(self, event):
        self.running = not self.running

    def _on_reset(self, event):
        self.running = False
        self._initialize_simulation()
        self._update_display()
        self.fig.canvas.draw_idle()

    def _on_step(self, event):
        if self.sim:
            self.sim.step()
            self._refresh_colors()
            self._update_display()
            self.fig.canvas.draw_idle()

    def _on_bound_toggle(self, label):
        self.params.bounding_enabled = not self.params.bounding_enabled
        if self.sim:
            self.sim.params.bounding_enabled = self.params.bounding_enabled

    def _initialize_simulation(self):
        """Initialize or reset the simulation."""
        self.sim = CentromereSimulator(self.params.copy(), seed=None)
        self.sim.initialize()
        self._refresh_colors()

    def _refresh_colors(self):
        """Refresh the colorizer with current unique sequences."""
        if self.sim:
            unique_seqs = self.sim.get_unique_sequences()
            self.colorizer.clear_cache()
            self.colorizer.colorize_batch(unique_seqs)

    def _animation_step(self, frame):
        """Animation callback."""
        if self.running and self.sim:
            self.sim.step()

            # Refresh colors periodically (every 10 generations)
            if self.sim.state.generation % 10 == 0:
                self._refresh_colors()

            self._update_display()

        return [self.image, self.stats_text]

    def show(self):
        """Display the visualization."""
        self._setup_figure()
        self._initialize_simulation()
        self._update_display()

        # Create animation
        self.anim = FuncAnimation(
            self.fig,
            self._animation_step,
            interval=int(self.generation_delay * 1000),
            blit=False,
            cache_frame_data=False
        )

        plt.show()


def main():
    """Run the visualizer."""
    print("Starting Centromere Evolution Visualizer...")
    print("Controls:")
    print("  - Start/Stop: Toggle simulation")
    print("  - Step: Advance one generation")
    print("  - Reset: Restart simulation")
    print("  - Sliders: Adjust mutation rates")
    print("  - Bounding: Toggle array size limits")
    print()

    viz = SimulationVisualizer(grid_width=100)
    viz.show()


if __name__ == "__main__":
    main()

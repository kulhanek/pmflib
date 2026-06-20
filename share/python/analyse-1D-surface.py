#!/usr/bin/env python3

import argparse
import math
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from scipy import ndimage
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize

# ------------------------------------------------------------------------------

rfac = 0.001987204258640       # kcal/mol/K

# ==============================================================================
# CV/Axis
# ==============================================================================

class Axis:
    """One collective-variable axis, represented internally on [0, 1]."""

    def __init__(self, cvmin, cvmax,nrbfs, npts, label):
        if cvmax <= cvmin:
            raise ValueError("cvmax must be larger than cvmin")
        if nrbfs <= 0:
            raise ValueError("nrbfs must be positive")
        if npts <= 0:
            raise ValueError("npts must be positive")

        self.cvmin = float(cvmin)
        self.cvmax = float(cvmax)
        self.nrbfs = int(nrbfs)
        self.npts = int(npts)
        self.label = label

        self.range = self.cvmax - self.cvmin
        self.width = 1.0 / float(self.nrbfs)
        self.width_scale = 1.0

        # Non-periodic CV: centers include both boundaries.
        self.centers = np.arange(self.nrbfs + 1, dtype=float) / float(self.nrbfs)

    # -------------------------------------------------------------------------

    def scale(self, x):
        """Convert a physical CV value to a scaled coordinate."""
        u = (np.asarray(x, dtype=float) - self.cvmin) / self.range
        return u

    # -------------------------------------------------------------------------

    def unscale(self, u):
        """Convert a scaled coordinate to a physical CV value."""
        u = np.asarray(u, dtype=float)
        return self.cvmin + u * self.range

    # -------------------------------------------------------------------------

    def delta(self, u):
        """
        Difference between scaled coordinate u and all RBF centers.
        """
        du = float(u) - self.centers
        return du

# ==============================================================================
# EnergySurface1D
# ==============================================================================

class EnergySurface1D:
    """
    One-dimensional energy surface represented by Gaussian radial basis functions.

    Internal coordinate:

        u = (x - xmin) / (xmax - xmin)                                      (1)

    RBF expansion:

        E(u) = sum_i A_i exp[-1/2 ((u - u_i) / s)^2]                         (2)

    The method eval(x) returns derivatives with respect to the original
    coordinate x. The method eval_u(u) returns derivatives with respect to the
    scaled coordinate u.
    """

    def __init__(
        self,
        cv1min,
        cv1max,
        cv1_nrbfs,
        cv1_nbins,
        cv1_label,
        ene_label,
        zmax,
        thrfac,
        temp,
        random_seed,
    ):
        self.x_axis = Axis(cv1min, cv1max, cv1_nrbfs, cv1_nbins, cv1_label)

        self.x_data = None
        self.e_data = None
        self.e_fit = None
        self.fit_residuals = None
        self.fit_rmse = None

        self.x_grid = None
        self.u_grid = None
        self.ENE = None
        self.SNG = None
        self.unsampled_mask = None
        self.sampled_regions = []
        self.sampled_region_map = None

        self.temp = float(temp)
        self.thr = float(thrfac) * self.temp * rfac

        self.ene_label = ene_label
        self.zmax = float(zmax)

        self.amplitudes = None
        self.sp_guesses = []
        self.sp_optimized = []
        self.basins = []

        if random_seed is None:
            random_seed = int(time.time())
        print(f"  Random seed: {random_seed}")
        self.rng = np.random.default_rng(random_seed)

    # --------------------------------------------------------------------------
    # Coordinate conversion
    # --------------------------------------------------------------------------

    def to_scaled(self, x):
        """Convert a physical coordinate to a scaled coordinate."""
        return float(self.x_axis.scale(float(x)))

    # --------------------------------------------------------------------------

    def from_scaled(self, u):
        """Convert a scaled coordinate to a physical coordinate."""
        return float(self.x_axis.unscale(float(u)))

    # --------------------------------------------------------------------------
    # Data loading and fitting
    # --------------------------------------------------------------------------

    def load(self, filename, xcolumn=1, ecolumn=2):
        """
        Load input data from a text file.

        Comments and empty lines are ignored.
            xcolumn: cv1
            ecolumn: energy
        """
        x_values = []
        e_values = []

        with open(filename, "r", encoding="utf-8") as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                line = line.split("#", 1)[0].strip()
                if not line:
                    continue

                fields = line.split()
                if len(fields) < max(xcolumn,ecolumn):
                    continue

                x_values.append(float(fields[xcolumn-1]))
                e_values.append(float(fields[ecolumn-1]))

        if len(e_values) == 0:
            raise ValueError("No valid data points were loaded")

        self.x_data = np.asarray(x_values, dtype=float)
        self.e_data = np.asarray(e_values, dtype=float)

        emin = np.nanmin(self.e_data)
        self.e_data = self.e_data - emin

    # --------------------------------------------------------------------------

    def _basis_1d(self, u):
        """
        Return 1D Gaussian basis values and derivatives in scaled coordinates.

        Returns
        -------
        g : ndarray
            Gaussian values.
        dg : ndarray
            First derivatives with respect to u.
        d2g : ndarray
            Second derivatives with respect to u.
        """
        du = self.x_axis.delta(u)
        s = self.x_axis.width * self.x_axis.width_scale

        g = np.exp(-0.5 * (du / s) ** 2)
        dg = g * (-du / s**2)
        d2g = g * ((du**2 / s**4) - (1.0 / s**2))

        return g, dg, d2g

    # --------------------------------------------------------------------------

    def _build_design_matrix(self, indices):
        """
        Build the RBF design matrix for selected data points.
        """
        
        indices = np.asarray(indices, dtype=int)
        B = np.empty((len(indices), len(self.x_axis.centers)), dtype=float)

        for row, k in enumerate(indices):
            u = self.to_scaled(self.x_data[k])
            B[row, :] = self._basis_1d(u)[0]

        return B

    # --------------------------------------------------------------------------

    @staticmethod
    def _solve_svd(B, y, rcond):
        """
        Solve a linear least-squares problem using an SVD pseudoinverse.
        """

        U, s, Vt = np.linalg.svd(B, full_matrices=False)
        if len(s) == 0:
            raise RuntimeError("SVD failed: no singular values found")

        cutoff = rcond * np.max(s)
        sinv = np.zeros_like(s)
        mask = s > cutoff
        sinv[mask] = 1.0 / s[mask]

        return Vt.T @ (sinv * (U.T @ y))

    # --------------------------------------------------------------------------

    def fit(self, sx=1.0, rcond=1.0e-12):
        """
        Fit RBF amplitudes to loaded data using an SVD pseudoinverse.
        """
        
        if self.x_data is None or self.e_data is None:
            raise RuntimeError("No data loaded. Call load() first.")

        self.x_axis.width_scale = float(sx)

        B = self._build_design_matrix(np.arange(len(self.e_data)))
        coeff = self._solve_svd(B, self.e_data, rcond)
        self.amplitudes = coeff.copy()

        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

        return self.amplitudes

    # --------------------------------------------------------------------------

    def fit_width_optimize(
        self,
        rcond=1.0e-12,
        width_scales=(0.5, 2.5),
        validation_fraction=0.2,
        verbose=False,
    ):
        """
        Fit RBF amplitudes and optimise the Gaussian width by scalar minimisation.
        """
        
        if self.x_data is None or self.e_data is None:
            raise RuntimeError("No data loaded. Call load() first.")

        if rcond <= 0.0:
            raise ValueError("rcond must be positive")
        if not (0.0 < validation_fraction < 1.0):
            raise ValueError("validation_fraction must be between 0 and 1")

        sx_min, sx_max = map(float, width_scales)
        if sx_min <= 0.0:
            raise ValueError("width-scale bounds must be positive")
        if sx_max <= sx_min:
            raise ValueError("width_scales must be given as (lower, upper)")

        ndata = len(self.e_data)
        indices = np.arange(ndata)
        self.rng.shuffle(indices)

        nvalid = max(1, int(round(validation_fraction * ndata)))
        valid_idx = indices[:nvalid]
        train_idx = indices[nvalid:]
        if len(train_idx) == 0:
            raise ValueError("not enough data points for train/validation split")

        def objective(sx):
            if sx <= 0.0:
                return np.inf

            self.x_axis.width_scale = float(sx)
            try:
                B_train = self._build_design_matrix(train_idx)
                coeff = self._solve_svd(B_train, self.e_data[train_idx], rcond)
                B_valid = self._build_design_matrix(valid_idx)
                e_pred = B_valid @ coeff
                rmse = float(np.sqrt(np.mean((e_pred - self.e_data[valid_idx]) ** 2)))
            except Exception:
                rmse = np.inf

            if verbose:
                print(f"  sx = {sx:12.6f}  RMSE = {rmse:12.6f}")

            return rmse

        result = minimize_scalar(objective, method="bounded", bounds=(sx_min, sx_max))
        if not result.success and verbose:
            print("Width optimisation warning:", result.message)

        self.x_axis.width_scale = float(result.x)
        B = self._build_design_matrix(np.arange(ndata))
        coeff = self._solve_svd(B, self.e_data, rcond)
        self.amplitudes = coeff.copy()

        self.best_width_scale_x = float(result.x)
        self.best_validation_rmse = float(result.fun)
        self.width_optimization_result = result
        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

        return self.amplitudes

    # --------------------------------------------------------------------------
    # Evaluation
    # --------------------------------------------------------------------------

    def eval_u(self, u):
        """
        Evaluate the surface in scaled coordinates.

        Returns value, first derivative and second derivative with respect to u.
        """
        if self.amplitudes is None:
            raise RuntimeError("The surface has not been fitted. Call fit() first.")

        g, dg, d2g = self._basis_1d(float(u))
        value = float(np.dot(self.amplitudes, g))
        gradient = float(np.dot(self.amplitudes, dg))
        hessian = float(np.dot(self.amplitudes, d2g))

        return value, gradient, hessian

    # --------------------------------------------------------------------------

    def eval(self, x):
        """
        Evaluate the surface in original coordinates.

        Returns value, first derivative and second derivative with respect to x.
        """
        u = self.to_scaled(x)
        value, gradient_u, hessian_u = self.eval_u(u)
        rx = self.x_axis.range
        gradient = gradient_u / rx
        hessian = hessian_u / (rx * rx)

        return value, gradient, hessian

    # --------------------------------------------------------------------------

    def eval_u_sng_f(self, u):
        """
        Stationary-point objective: squared norm of energy gradient in scaled coordinates.
        """
        
        value, gradient, hessian = self.eval_u(float(u))
        return gradient * gradient
    
    # --------------------------------------------------------------------------

    def eval_u_sng_fg(self, u):
        """
        Stationary-point objective: squared norm of energy gradient in scaled coordinates.
        """
        
        value, gradient, hessian = self.eval_u(float(u))

        sng_f = gradient * gradient
        sng_g = 2.0 * gradient * hessian

        return sng_f, np.array([sng_g],dtype=float)
    
    # --------------------------------------------------------------------------

    def calc_ene_and_sng(self):
        """
        Evaluate the ENE profile and SNG metric on a regular 1D grid.
        """

        if self.x_axis.npts <= 1:
            raise ValueError("cv1nbins must be larger than 1")

        x_edges = np.linspace(self.x_axis.cvmin, self.x_axis.cvmax, self.x_axis.npts + 1)
        self.x_grid = 0.5 * (x_edges[:-1] + x_edges[1:])

        u_edges = np.linspace(0.0, 1.0, self.x_axis.npts + 1)
        self.u_grid = 0.5 * (u_edges[:-1] + u_edges[1:])

        self.ENE = np.empty_like(self.x_grid, dtype=float)
        self.SNG = np.empty_like(self.x_grid, dtype=float)

        for i, x in enumerate(self.x_grid):
            self.ENE[i], _, _ = self.eval(x)
            self.SNG[i] = self.eval_u_sng_f(self.u_grid[i])

        self.unsampled_mask, nearest_dist = self.detect_unsampled_grid_points()
        self.ENE[self.unsampled_mask] = self.zmax

        spmax = np.nanmax(self.SNG)
        self.SNG[self.unsampled_mask] = spmax

    # --------------------------------------------------------------------------

    def _default_sample_spacing_u(self):
        """
        Return a robust nearest-neighbour spacing estimate in scaled units.
        """

        if self.x_data is None or len(self.x_data) < 2:
            return 1.0 / float(self.x_axis.npts)

        data_u = np.sort(np.asarray(self.x_axis.scale(self.x_data), dtype=float))
        data_u = np.unique(data_u)

        if len(data_u) < 2:
            return 1.0 / float(self.x_axis.npts)

        diffs = np.diff(data_u)

        diffs = diffs[diffs > 0.0]
        if len(diffs) == 0:
            return 1.0 / float(self.x_axis.npts)

        return float(np.min(diffs))

    # --------------------------------------------------------------------------

    def build_sampled_regions(self, max_distance):
        """
        Build connected sampled regions from the real sampled coordinates.

        The grid is used only to decide whether a grid point is sampled or
        unsampled.  Basin limits are taken from the actual sampled data points
        stored in ``self.x_data``.  Each sampled grid point is mapped to one
        dictionary item with at least ``min`` and ``max`` physical sampled
        coordinates.
        """

        data_x = np.asarray(self.x_data, dtype=float)
        data_u = np.asarray(self.x_axis.scale(data_x), dtype=float)

        if len(data_u) == 0:
            self.sampled_regions = []
            self.sampled_region_map = np.empty_like(self.x_grid, dtype=object)
            self.sampled_region_map[:] = None
            return

        order = np.argsort(data_u)
        sorted_u = data_u[order]
        sorted_x = data_x[order]

        # Two neighbouring sampled points belong to different sampled regions
        # only if the gap between them is wide enough to contain grid points
        # whose nearest sampled point is farther than max_distance.
        split_distance = 2.0 * float(max_distance)

        regions = []
        sorted_index_to_region = {}

        def add_region(component):
            rid = len(regions) + 1
            comp = np.asarray(component, dtype=int)
            comp_u = sorted_u[comp]
            comp_x = sorted_x[comp]

            item = {
                "id": rid,
                "min": float(np.min(comp_x)),
                "max": float(np.max(comp_x)),
                "umin": float(np.min(comp_u)),
                "umax": float(np.max(comp_u)),
            }
            regions.append(item)
            for idx in comp:
                sorted_index_to_region[int(idx)] = item

        start = 0
        for i, gap in enumerate(np.diff(sorted_u)):
            if gap > split_distance:
                add_region(np.arange(start, i + 1))
                start = i + 1
        add_region(np.arange(start, len(sorted_u)))
        
        grid_u = np.asarray(self.x_axis.scale(self.x_grid), dtype=float)
        region_map = np.empty(len(grid_u), dtype=object)
        region_map[:] = None

        for i, gu in enumerate(grid_u):
            nearest_sorted_idx = int(np.argmin(np.abs(gu - sorted_u)))
            nearest_dist = float(abs(gu - sorted_u[nearest_sorted_idx]))

            if nearest_dist <= max_distance:
                region_map[i] = sorted_index_to_region[nearest_sorted_idx]

        self.sampled_regions = regions
        self.sampled_region_map = region_map

    # --------------------------------------------------------------------------

    def detect_unsampled_grid_points(self, max_distance=None):
        """
        Detect grid points too far from any sampled data point.
        """
        
        if self.x_data is None:
            raise RuntimeError("No data loaded. Call load() first.")
        
        if self.x_grid is None:
            raise RuntimeError("Regular grid is not available. Call calc_ene_and_sng() first.")

        data_u = np.asarray(self.x_axis.scale(self.x_data), dtype=float).reshape(-1, 1)
        grid_u = np.asarray(self.x_axis.scale(self.x_grid), dtype=float).reshape(-1, 1)

        if max_distance is None:
            # max_distance is a distance in scaled coordinates.
            max_distance = 1.05 * self._default_sample_spacing_u()

        # Direct computation is sufficient and avoids an extra dependency.
        nearest_dist = np.min(np.abs(grid_u - data_u.T), axis=1)

        unsampled_mask = nearest_dist > max_distance
        self.build_sampled_regions(max_distance)
        return unsampled_mask, nearest_dist

    # --------------------------------------------------------------------------
    # Plot helpers
    # --------------------------------------------------------------------------

    @staticmethod
    def make_gnuplot_like_colormaps():
        """
        Create colormaps compatible with the common PMFLib/gnuplot palette.
        """

        gnuplot_colors = [
            "#3288BD",
            "#66C2A5",
            "#ABDDA4",
            "#E6F598",
            "#FEE08B",
            "#FDAE61",
            "#F46D43",
            "#D53E4F",
        ]

        basin_cmap = ListedColormap(gnuplot_colors, name="gnuplot_like_regions")

        return basin_cmap
    
    # --------------------------------------------------------------------------

    def _finish_plot(self, fig, ax, show, save, dpi):

        ax.axhline(0.0, linestyle="-", linewidth=0.8, zorder=0, color="black")

        ax.set_xlabel(self.x_axis.label)
        ax.set_ylabel(self.ene_label)
        ax.set_xlim(self.x_axis.cvmin, self.x_axis.cvmax)
        ax.set_ylim(bottom=-1.0)
        fig.tight_layout()

        if save is not None:
            fig.savefig(save, dpi=dpi)
        if show:
            plt.show()

        plt.close()

    # --------------------------------------------------------------------------

    def plot_raw(self, marker_size=25, show=True, save=None, figsize=None, dpi=300):
        """Plot loaded scattered 1D energy data."""

        if self.e_data is None:
            raise RuntimeError("The surface has not been loaded. Call load() first.")

        fig, ax = plt.subplots(figsize=figsize)

        ax.axhline(self.zmax, linestyle="--", zorder=0, linewidth=0.8)
        ax.scatter(self.x_data, self.e_data, marker="x", zorder=1, s=6, alpha=0.7, color="red", label="data")

        ax.set_title("Original 1D FES")
        self._finish_plot(fig, ax, show, save, dpi)
    
    # --------------------------------------------------------------------------

    def plot_rbf(self, show=True, save=None, figsize=None, dpi=300):
        """Plot the fitted 1D energy profile."""

        if self.ENE is None:
            raise RuntimeError("The interpolated profile is not available. Call calc_ene_and_sng() first.")

        fig, ax = plt.subplots(figsize=figsize)
        
        ax.plot(self.x_grid, np.minimum(self.ENE, self.zmax), zorder=0, linewidth=1.0, label="RBF")
        ax.axhline(self.zmax, linestyle="--", zorder=0, linewidth=0.8)

        if self.x_data is not None:
            ax.scatter(self.x_data, np.minimum(self.e_data, self.zmax), marker="x", zorder=1, s=6, alpha=0.7, color="red", label="data")

        ax.legend(loc="best")

        ax.set_title("Interpolated 1D FES")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_u_sng(self, zmax=None, show=True, save=None, figsize=None, dpi=300):
        """Plot the stationary-point metric on the scaled coordinate."""

        if self.SNG is None:
            raise RuntimeError("The SNG metric is not available. Call calc_ene_and_sng() first.")

        fig, ax = plt.subplots(figsize=figsize)

        sp_plot = self.SNG if zmax is None else np.minimum(self.SNG, float(zmax))

        ax.plot(self.u_grid, sp_plot, linewidth=1.0)
        ax.set_xlabel("CV1 [scaled]")
        ax.set_ylabel("SP metric")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(bottom=0.0)
        ax.set_title("SP Metric")
        fig.tight_layout()

        if save is not None:
            fig.savefig(save, dpi=dpi)
        if show:
            plt.show()

        plt.close()

    # --------------------------------------------------------------------------

    def plot_rbf_with_sp_guesses(self, show=True, save=None, figsize=None, dpi=300):
        """Plot the fitted 1D profile with stationary-point guesses."""

        if self.ENE is None:
            raise RuntimeError("The interpolated profile is not available. Call calc_profile_and_sp() first.")

        fig, ax = plt.subplots(figsize=figsize)

        ax.plot(self.x_grid, np.minimum(self.ENE, self.zmax), linewidth=1.0, label="RBF",zorder=0)
        ax.axhline(self.zmax, linestyle="--", linewidth=0.8,zorder=0)

        if len(self.sp_guesses) > 0:
            xs = [p["x"] for p in self.sp_guesses]
            ys = [self.eval(p["x"])[0] for p in self.sp_guesses]
            ax.scatter(xs, ys, s=35, marker="o", facecolors="none", edgecolors="orange", linewidths=1.0,zorder=1)
            ax.scatter(xs, ys, s=10, marker="x", color="orange",zorder=1)
            for p, y in zip(self.sp_guesses, ys):
                ax.text(p["x"], y, f" {p['id']}", color="orange", fontsize=9, va="center",zorder=2)

        ax.set_title("Interpolated 1D FES with Stationary-Point Guesses")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_rbf_with_sp_optimized(self, show=True, save=None, figsize=None, dpi=300):
        """Plot the fitted 1D profile with optimized stationary points."""

        if self.ENE is None:
            raise RuntimeError("The interpolated profile is not available. Call calc_profile_and_sp() first.")

        fig, ax = plt.subplots(figsize=figsize)
        
        ax.plot(self.x_grid, np.minimum(self.ENE, self.zmax), linewidth=1.0, label="RBF",zorder=0)
        ax.axhline(self.zmax, linestyle="--", linewidth=0.8,zorder=0)

        colors = {
                    "T": "red",
                    "S": "blue",
                }

        for p in self.sp_optimized:
            x = p["x"]
            y = p["ene"]
            ax.scatter([x], [y], s=35, marker="o", facecolors="none", edgecolors=colors[p["type"]], linewidths=1.0,zorder=1)
            ax.scatter([x], [y], s=10, marker="x", color=colors[p["type"]], linewidths=1.0,zorder=1)
            ax.text(x, y, f" {p['id']}:{p['type']}", color=colors[p["type"]], fontsize=9, va="center",zorder=2)

        ax.set_title("Interpolated 1D FES with Optimized Stationary Points")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_rbf_with_basins(self, show=True, save=None, figsize=None, dpi=300):
        """Plot the fitted 1D profile with minimum basins as shaded intervals."""

        if self.ENE is None:
            raise RuntimeError("The interpolated profile is not available. Call calc_profile_and_sp() first.")

        ws_cmap = self.make_gnuplot_like_colormaps()

        fig, ax = plt.subplots(figsize=figsize)
        
        ax.plot(self.x_grid, np.minimum(self.ENE, self.zmax), linewidth=1.0, label="RBF",zorder=0)
        ax.axhline(self.zmax, linestyle="--", linewidth=0.8,zorder=0)

        colors = {
                    "T": "red",
                    "S": "blue",
                }

        for p in self.sp_optimized:
            x = p["x"]
            y = p["ene"]
            ax.scatter([x], [y], s=35, marker="o", facecolors="none", edgecolors=colors[p["type"]], linewidths=1.0,zorder=1)
            ax.scatter([x], [y], s=10, marker="x", color=colors[p["type"]], linewidths=1.0,zorder=1)
            ax.text(x, y, f" {p['id']}:{p['type']}", color=colors[p["type"]], fontsize=9, va="center", ha="left",zorder=2)

        if len(self.basins) > 0:
            ymin, ymax = ax.get_ylim()
            for bidx, basin in enumerate(self.basins):
                color = ws_cmap(bidx % ws_cmap.N)
                xmin = basin["xmin"]
                xmax = basin["xmax"]
                ax.axvspan(xmin, xmax, alpha=0.25, color=color, zorder=-1)

                ax.text(
                    0.5 * (xmin + xmax),
                    0.95 * ymax,
                    str(basin["id"]),
                    ha="center",
                    va="top",
                    fontsize=9,
                    zorder=0
                )
            ax.set_ylim(ymin, ymax)

        ax.set_title("Interpolated 1D FES with Minima Basins")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------
    # Stationary-point detection and optimisation
    # --------------------------------------------------------------------------

    def detect_sp_regions_connected(self, eps, min_size=1, representative="minabs"):
        """
        Detect connected 1D intervals where the squared derivative is small.
        """
        
        if self.SNG is None or self.ENE is None:
            raise RuntimeError("The SNG metric is not available. Call calc_ene_and_sng() first.")

        self.sp_guess_mask = (self.SNG >= 0.0) & (self.SNG <= eps) & (self.ENE < self.zmax)
        self.sp_guess_labels, nlabels = ndimage.label(self.sp_guess_mask)

        self.sp_guesses = []
        for label in range(1, nlabels + 1):
            idx = np.argwhere(self.sp_guess_labels == label).ravel()
            if len(idx) < min_size:
                self.sp_guess_labels[self.sp_guess_labels == label] = 0
                continue

            values = self.SNG[idx]
            if representative == "minabs":
                k = int(idx[np.argmin(np.abs(values))])
            elif representative == "centroid":
                centroid = idx.mean()
                k = int(idx[np.argmin((idx - centroid) ** 2)])
            else:
                raise ValueError("representative must be 'minabs' or 'centroid'")

            self.sp_guesses.append({
                "id": int(label),
                "i": int(k),
                "u": float(self.u_grid[k]),
                "sng": float(self.SNG[k]),
                "x": float(self.x_grid[k]),
                "ene": float(self.ENE[k]),
                "size": int(len(idx)),
            })

        return self.sp_guesses, self.sp_guess_labels, self.sp_guess_mask

    # --------------------------------------------------------------------------

    def find_sp(self, gx, trust_r0, rsize = None):
        """
        Locate and classify a 1D stationary point near the initial guess.
        """

        # Input guess
        u0s = np.array([float(self.to_scaled(gx))], dtype=float)

        # -------------------------------------------------------------------------
        # Find SP
        # -------------------------------------------------------------------------

        best_res = None
        best_fun = np.inf
        best_u0 = None

        print("")

        # we need to search in a small box in the case of very narrow minima
        if rsize is not None:
            # used trust region detected by detect_sp_regions_connected
            trust_r = float(rsize) * 1.0 / self.x_axis.npts
        else:
            trust_r = trust_r0

        bounds = [(u0s[0]-trust_r, u0s[0]+trust_r)]

        print("Stationary point scan (L-BFGS-B):")
        # try gradient based optimizer
        for attempt in range(0, 10):

            if attempt > 0:
                dx = self.rng.uniform(
                    low=-trust_r/3.0,  # use smaller interval
                    high=trust_r/3.0,
                    size=u0s.shape
                )

                u0 = np.zeros(1)
                u0 = u0s + dx
            else:
                u0 = u0s

            print(f"  u0 = {u0}")

            res = minimize(
                self.eval_u_sng_fg,
                u0,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds
            )

            print(f"      err = {res.fun:12.6e}, uopt = {res.x}")

            if (res.x[0] < 0.0) or (res.x[0] > 1.0):
                print(" >> solution out-of-box: ignoring")
                continue

            if np.linalg.norm(res.x-u0) > trust_r:
                print(" >> solution out-of-trust-region: ignoring")
                continue 

            if res.fun < best_fun:
                best_fun = res.fun
                best_res = res
                best_u0 = u0.copy()

        if best_res is None:
            print("Stationary point scan Nelder-Mead):")
            # try non-gradient optimizers
            for attempt in range(0, 10):

                if attempt > 0:
                    du = self.rng.uniform(
                        low=-trust_r/3.0,  # use smaller interval
                        high=trust_r/3.0,
                        size=u0s.shape
                    )

                    u0 = np.zeros(2)
                    u0 = u0s + du
                else:
                    x0 = u0s

                print(f"  x0 = {x0}")

                res = minimize(
                    self.eval_u_sng_f,
                    u0,
                    method="Nelder-Mead",
                    bounds=bounds
                )

                print(f"      err = {res.fun:12.6e}, uopt = {res.x}")

                if (res.x[0] < 0.0) or (res.x[0] > 1.0):
                    print(" >> solution out-of-box: ignoring")
                    continue

                if np.linalg.norm(res.x-u0) > trust_r:
                    print(" >> solution out-of-trust-region: ignoring")
                    continue 

                if res.fun < best_fun:
                    best_fun = res.fun
                    best_res = res
                    best_u0 = u0.copy()


        if best_res is None:
            # not found
            item = {
                'id': -1
            }
            return item

        print("")
        print("Optimized SP parameters:")
        print(f"  initial x0      = {best_u0}")
        print(f"  best error      = {best_fun:12.6e}")
        print(f"  best parameters = {best_res.x}")

        print("")
        print(best_res)

        uopt = best_res.x[0]             # scaled
        xopt = self.from_scaled(uopt) 

        print("")

        ene, _, hess_u = self.eval_u(uopt)
        sng = self.eval_u_sng_f(uopt)

        if hess_u > 0.0:
            sp_type = "S"
            sp_type_text = "local minimum"
        elif hess_u < 0.0:
            # In a 1D free-energy profile this is the analogue of a transition state.
            sp_type = "T"
            sp_type_text = "local maximum / transition state"
        else:
            sp_type = "F"
            sp_type_text = "flat/inconclusive"

        print(f"  uopt = {uopt:12.6f}")
        print(f"  xopt = {xopt:12.6f}")
        print(f"  ene  = {ene:12.6f}")
        print(f"  sng  = {sng:12.6e}")
        print(f"  sl1  = {hess_u:12.6f}")
        print(f"  Type of SP: {sp_type_text}")

        # Characteristic half-width of the quadratic neighbourhood at self.thr.
        if abs(hess_u) > 0.0:
            r_scaled = math.sqrt(max(0.0, 2.0 * self.thr / abs(hess_u)))
            r_original = r_scaled * self.x_axis.range
        else:
            r_original = np.inf

        item = {
            "id": len(self.sp_optimized) + 1,
            "x": float(xopt),
            "ene": float(ene),
            "sng": float(sng),
            "r1": float(r_original),
            "type": sp_type,
            "sl1": float(hess_u),
        }

        self.sp_optimized.append(item)
        return item

    # --------------------------------------------------------------------------

    def remove_near_duplicate_sp_optimized(self, min_distance_u=0.01):
        """
        Remove nearly duplicate optimized stationary points.

        Two points are considered duplicates only if they have the same type
        ('S', 'T', or 'M') and their distance in scaled UV coordinates is
        smaller than or equal to min_distance_u.

        From each duplicate group, the point with the lowest energy is kept.

        Parameters
        ----------
        min_distance_u : float
            Minimum distance between two SPs to be considered as individual points.

        Returns
        -------
        removed : list of dict
            List of removed stationary points.
        """

        print("")
        print("# Detecting near duplicate points ...")

        if min_distance_u <= 0.0:
            raise ValueError("min_distance_u must be positive")

        kept = []
        removed = []

        def better(a, b):
            if a["type"] == "S":
                return a["ene"] <= b["ene"]
            if a["type"] == "T":
                return a["ene"] >= b["ene"]
            return True

        for p in sorted(self.sp_optimized, key=lambda item: item["x"]):
            p_type = p.get("type", "")
            if p_type not in ("S", "T"):
                kept.append(p)
                continue

            p_u = self.to_scaled(p["x"])
            duplicate_index = None

            for i, q in enumerate(kept):
                if q.get("type", "") != p_type:
                    continue

                q_u = self.to_scaled(q["x"])
                du = p_u - q_u

                if abs(du) <= min_distance_u:
                    duplicate_index = i
                    break

            if duplicate_index is None:
                kept.append(p)
            else:
                q = kept[duplicate_index]
                if better(p, q):
                    removed.append(q)
                    kept[duplicate_index] = p
                else:
                    removed.append(p)

        self.sp_optimized = sorted(kept, key=lambda item: int(item["id"]))

        print(f"  Removed {len(removed)} near-duplicate optimized stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return removed

    # --------------------------------------------------------------------------

    def remove_stationary_point_outliers(self, minslam=5.0, maxslam=5000.0):
        """Remove stationary points with implausible scaled Hessian values."""

        print("")
        print("# Detecting Hessian outliers ...")
        print(f"  Allowed Hessian eigenvalues: <{minslam:.1f}, {maxslam:.1f}>")

        kept = []
        removed = []

        for item in self.sp_optimized:
            outlier = abs(item["sl1"]) < minslam or abs(item["sl1"]) > maxslam
            if outlier:
                print(
                    f"  Removing stationary point {item['id']} ({item['type']}) "
                    f"due to Hessian value {item['sl1']:.1f}."
                )
                removed.append(item)
            else:
                kept.append(item)

        self.sp_optimized = kept
        print(f"  Removed {len(removed)} stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return removed
    
    # --------------------------------------------------------------------------

    def remove_sng_outliers(self, max_sng=0.5):

        print(f"")
        print(f"# Detecting SNG outliers ...")
        print(f"  Max sng: {max_sng:.2f}")
    
        kept = []
        removed = []

        for item in self.sp_optimized:
            outlier = (
                item["sng"] > max_sng
            )
            if outlier:
                print(f"  Removing stationary point {item['id']} ({item['type']}) due to large sng: ({item['sng']:.2f}) out-of-allowed values.")
                removed.append(item)
            else:
                kept.append(item)

        self.sp_optimized = kept

        print(f"  Removed {len(removed)} stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return removed

    # --------------------------------------------------------------------------
    # Basins
    # --------------------------------------------------------------------------

    def _sampled_regions_u(self):
        """Return connected sampled regions in scaled coordinates."""

        if self.u_grid is None or self.ENE is None:
            raise RuntimeError("The regular grid is not available. Call calc_ene_and_sng() first.")

        if self.sampled_region_map is None:
            self.build_sampled_regions(1.05 * self._default_sample_spacing_u())

        used_ids = set()
        regions = []
        for item in self.sampled_region_map:
            if item is None or item["id"] in used_ids:
                continue
            used_ids.add(item["id"])
            regions.append(item)

        return sorted(regions, key=lambda item: item["umin"])

    # --------------------------------------------------------------------------

    @staticmethod
    def _region_contains_u(u, region):
        """Return True if scaled coordinate u belongs to a sampled region."""

        return float(region["umin"]) <= float(u) <= float(region["umax"])

    # --------------------------------------------------------------------------

    def calculate_basins(self):
        """
        Calculate 1D minima basins.

        Each minimum basin is bounded by neighbouring transition states/maxima
        within the same connected sampled region.  If a maximum is missing on
        one side, the basin is clipped at the sampled-region boundary.
        """

        minima = sorted(
            [p for p in self.sp_optimized if p["type"] == "S"],
            key=lambda p: self.to_scaled(p["x"]),
        )
        maxima = sorted(
            [p for p in self.sp_optimized if p["type"] == "T"],
            key=lambda p: self.to_scaled(p["x"]),
        )

        if len(minima) == 0:
            print("  WARNING: No minima found; basins cannot be calculated.")
            self.basins = []
            return self.basins

        sampled_regions = self._sampled_regions_u()
        if len(sampled_regions) == 0:
            print("  WARNING: No sampled region found; basins cannot be calculated.")
            self.basins = []
            return self.basins

        self.basins = []

        for pt in minima:
            pt_u = self.to_scaled(pt["x"])

            containing_region = None
            for region in sampled_regions:
                if self._region_contains_u(pt_u, region):
                    containing_region = region
                    break

            if containing_region is None:
                print(f"  WARNING: Minimum {pt['id']} is outside sampled regions; skipping basin.")
                continue

            u0 = float(containing_region["umin"])
            u1 = float(containing_region["umax"])

            maxima_in_region = [
                self.to_scaled(p["x"])
                for p in maxima
                if u0 < self.to_scaled(p["x"]) < u1
            ]

            left_maxima = [u for u in maxima_in_region if u < pt_u]
            right_maxima = [u for u in maxima_in_region if u > pt_u]

            umin = max(left_maxima) if left_maxima else u0
            umax = min(right_maxima) if right_maxima else u1
            xmin = self.from_scaled(umin)
            xmax = self.from_scaled(umax)

            pt["xmin"] = float(xmin)
            pt["xmax"] = float(xmax)
            pt["umin"] = float(umin)
            pt["umax"] = float(umax)

            self.basins.append({
                "id": pt["id"],
                "xmin": float(xmin),
                "xmax": float(xmax),
                "umin": float(umin),
                "umax": float(umax),
                "sampled_region": containing_region,
                "pt": pt,
            })

        return self.basins

    # --------------------------------------------------------------------------

    def _basin_grid_mask(self, basin):
        """Return a mask selecting sampled grid points belonging to one basin."""

        mask = (self.x_grid >= basin["xmin"]) & (self.x_grid <= basin["xmax"])

        if self.unsampled_mask is not None:
            mask &= ~self.unsampled_mask
        mask &= self.ENE < self.zmax
        return mask

    # --------------------------------------------------------------------------

    def calculate_basins_ene_state(self):
        """Calculate Boltzmann state energies for 1D basins."""

        minima = [pt for pt in self.sp_optimized if pt["type"] == "S"]

        if len(minima) == 0:
            return

        dx = self.x_axis.range / self.x_axis.npts

        for basin in self.basins:
            pt = basin["pt"]
            mask = self._basin_grid_mask(basin)
            if not np.any(mask):
                print(f"  WARNING: Basin {pt['id']} contains no sampled grid points; state energy is undefined.")
                pt["ene_state"] = float("nan")
                continue
            weights = np.exp(-self.ENE[mask] / (self.temp * rfac))
            Q = np.sum(weights) * dx
            ene_state = -self.temp * rfac * math.log(Q)
            pt["ene_state"] = float(ene_state)

        # get minimum energy and calculate corrected energy value
        valid_minima = [pt for pt in self.sp_optimized if pt['type'] == "S" and np.isfinite(pt.get('ene_state', np.nan))]
        if len(valid_minima) == 0:
            return

        ene_min = min(pt['ene'] for pt in valid_minima)
        ene_state_min = min(pt['ene_state'] for pt in valid_minima)
        for pt in self.sp_optimized:
            if pt['type'] == "S":
                pt['ene0'] = pt['ene'] - ene_min
                pt['ene_state0'] = pt['ene_state'] - ene_state_min if np.isfinite(pt.get('ene_state', np.nan)) else float("nan")

# ==============================================================================
# Command-line arguments
# ==============================================================================

def parse_args():

    def str2bool(value):
        """Parse robust command-line boolean values."""
        if isinstance(value, bool):
            return value

        value = value.lower()
        if value in ("yes", "true", "t", "1", "on"):
            return True
        if value in ("no", "false", "f", "0", "off"):
            return False

        raise argparse.ArgumentTypeError(f"Cannot interpret '{value}' as a boolean value.")

    def parse_figsize(value):
        """Converts a comma-separated string into a tuple of floats."""
        try:
            # Split the string by comma and convert to floats
            parts = value.split(',')
            if len(parts) != 2:
                raise ValueError()
            return tuple(map(float, parts))
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid figsize format: '{value}'. Must be 'width,height' (e.g., '10,6')."
            )
        
    # --------------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description="Detect and optimize stationary points on a 1D free-energy/potential-energy surface."
    )

    # --------------------------------------------------------------------------
    # CV1
    # --------------------------------------------------------------------------

    cv1group = parser.add_argument_group("The collective variable specification")

    cv1group.add_argument("--cv1label", type=str, default=r"cv1", 
        help="Label of the collective variable.")
    
    cv1group.add_argument("--cv1min", type=float, required=True, 
        help="Minimum value of the collective variable.")
    
    cv1group.add_argument("--cv1max", type=float, required=True, 
        help="Maximum value of the collective variable.")
    
    cv1group.add_argument("--cv1nbins", type=int, required=True, 
        help="Number of grid bins/points.")

    # --------------------------------------------------------------------------
    # Energy label
    # --------------------------------------------------------------------------

    enegroup = parser.add_argument_group("The energy axis specification")

    enegroup.add_argument("--enelabel", type=str, 
        default=r"${\Delta}G [kcal/mol]$", help="Energy label.")
    
    enegroup.add_argument("--zmax", type=float, required=True,
        help="Maximum energy value considered.")

    # --------------------------------------------------------------------------
    # RBF interpolation
    # --------------------------------------------------------------------------

    rbfgroup = parser.add_argument_group("The RBF (Radial Basis Function) interpolation specification")

    rbfgroup.add_argument("--cv1nrbfs", type=int, default=20, 
        help="Number of RBFs.")
    
    rbfgroup.add_argument("--rbfwidthmode", type=str, default="static", 
        help="RBF width mode: static, gridsearch, optimize.")
    
    rbfgroup.add_argument("--rbfsx", type=float, default=1.5, 
        help="Width factor in static width mode.")
    
    rbfgroup.add_argument("--rcond", type=float, default=1.0e-9, 
        help="SVD cutoff for RBF fitting.")

    # --------------------------------------------------------------------------

    filegroup = parser.add_argument_group("The input/output files specification")

    filegroup.add_argument("--input-fes", type=str, dest="fname_input_fes", required=True,
        help="Input FES filename." )
    
    filegroup.add_argument("--input-fes-x-column", type=int, default=1,
        help="Index of x-column in the input FES file." )
    
    filegroup.add_argument("--input-fes-e-column", type=int, default=2,
        help="Index of e-column in the input FES file." )

    filegroup.add_argument("--sp-guesses", type=str, dest="fname_sp_guesses", default="gpts.txt",
        help="Output/Input filename for initial stationary-point guesses." )

    filegroup.add_argument("--sp-optimized", type=str, dest="fname_sp_optimized", default="opts.txt",
        help="Output filename for optimized stationary points." )
    
    filegroup.add_argument("--sp-basins", type=str, dest="fname_sp_basins", default="bpts.txt",
        help="Output filename for basins stationary points." )
    
    filegroup.add_argument("--basins", type=str, dest="fname_basins", default="basins.txt",
        help="Output filename for basins." )
    
    # --------------------------------------------------------------------------

    sysgroup = parser.add_argument_group("The system specification")

    sysgroup.add_argument("--actions",type=str,default="all",
        help="Comma-separated actions: all, loadfes, optrbfs, calcsurfs, loadguess, guess, loadopts, optimize, basins.",
    )
    sysgroup.add_argument( "--temp", type=float, default=300.0,
        help="Thermodynamic temperature." )

    sysgroup.add_argument( "--thrfac", type=float, default=0.25,
        help="Energy threshold factor for defining ellipse representing a stationary point (thrfac * kB * T)." )

    sysgroup.add_argument( "--random_seed", type=int, default=None,
        help="Random generator seed." )


    # --------------------------------------------------------------------------
    # Thresholds
    # --------------------------------------------------------------------------

    thresholdgroup = parser.add_argument_group("Threshold specification")
    
    thresholdgroup.add_argument( "--sng-cutoff", type=float, default=20.0,
        help="Cut-off value of SNG to consider a stationary point." )

    thresholdgroup.add_argument( "--sng-min-size", type=int, default=1,
        help="Minimum number of bins to consider a SNG region as a stationary point." )
    
    thresholdgroup.add_argument( "--trust-r0", type=float, default=0.005,
        help="Trust radius in scaled coordinates to locate stationary points." )
    
    thresholdgroup.add_argument( "--minslam", type=float, default=5.0,
        help="Minimal allowed value for scaled Hessian eigenvalue." )

    thresholdgroup.add_argument( "--maxslam", type=float, default=50000.0,
        help="Maximum allowed value for scaled Hessian eigenvalue." )

    thresholdgroup.add_argument( "--min-distance-u", type=float, default=0.005,
        help="Minimum distance between two SPs to be considered as individual points." )

    thresholdgroup.add_argument( "--max-sng", type=float, default=0.5,
        help="Maximum value of SNG for optimized stationary points." )

    # --------------------------------------------------------------------------
    # Plots
    # --------------------------------------------------------------------------

    plotgroup = parser.add_argument_group("The graphical plot specification")

    plotgroup.add_argument("--showrawfes", type=str2bool, default=False,
        help="Show raw FES plot.")
    
    plotgroup.add_argument("--saverawfes", type=str, default="FigureFES-RAW.png",
        help="Save raw FES plot.")
    
    plotgroup.add_argument("--showrbffes", type=str2bool, default=False,
        help="Show RBF FES plot.")
    
    plotgroup.add_argument( "--showsng", type=str2bool, default=False,
        help="Show an interactive plot with the SNG regions (action: guess)." )

    plotgroup.add_argument( "--savesng", type=str, default="FigureFES-SNG.png",
        help="Save a plot with the SNG regions (action: guess)." )
    
    plotgroup.add_argument("--saverbffes", type=str, default="FigureFES-RBF.png",
        help="Save RBF FES plot.")
    
    plotgroup.add_argument("--showgpts", type=str2bool, default=False,
        help="Show FES with guessed SPs.")
    
    plotgroup.add_argument("--savegpts", type=str, default="FigureFES-GPTS.png",
        help="Save FES with guessed SPs.")
    
    plotgroup.add_argument("--showopts", type=str2bool, default=False,
        help="Show FES with optimized SPs.")
    
    plotgroup.add_argument("--saveopts", type=str, default="FigureFES-OPTS.png",
        help="Save FES with optimized SPs.")
    
    plotgroup.add_argument("--showbasins", type=str2bool, default=False,
        help="Show FES with minima basins.")
    
    plotgroup.add_argument("--savebasins", type=str, default="FigureFES-Basins.png",
        help="Save FES with minima basins.")
    
    plotgroup.add_argument('--figsize',type=parse_figsize,default=(6.4, 4.8),  # Default Matplotlib size fallback
        help="Figure size as 'width,height' in inches (default: 6.4,4.8)" )
    
    plotgroup.add_argument("--dpi", type=int, default=300,
        help="Resolution for plot figures.")

    return parser.parse_args()

# ==============================================================================
# Action handlers
# ==============================================================================

def load_fes(args, surf):
    print("")
    print(f"# Load FES: {args.fname_input_fes}")
    surf.load(args.fname_input_fes,xcolumn=args.input_fes_x_column,ecolumn=args.input_fes_e_column)

    if args.showrawfes or args.saverawfes is not None:
        surf.plot_raw(show=args.showrawfes, save=args.saverawfes, dpi=args.dpi)

# ------------------------------------------------------------------------------

def opt_rbfs(args, surf):
    print("")
    print("# Optimize RBF ...")
    print(f"  Width mode: {args.rbfwidthmode}")

    if args.rbfwidthmode == "static":
        print(f"  Sx:         {args.rbfsx:10.3f}")
        surf.fit(sx=args.rbfsx, rcond=args.rcond)
    elif args.rbfwidthmode == "optimize":
        surf.fit_width_optimize(rcond=args.rcond, verbose=True)
    else:
        raise ValueError(f"Unsupported RBF width mode: {args.rbfwidthmode:s}")

    print(f"  RMSE:       {surf.fit_rmse:10.3f}")

# ------------------------------------------------------------------------------

def cal_surfs(args, surf):
    print("")
    print("# Calculate ENE and SNG profiles ...")
    surf.calc_ene_and_sng()

    if args.showrbffes or args.saverbffes is not None:
        surf.plot_rbf(show=args.showrbffes, save=args.saverbffes, figsize=args.figsize, dpi=args.dpi)

# ------------------------------------------------------------------------------

def load_sp_guesses(args, surf):
    print("")
    print("# Load initial stationary point guesses (2nd column -> x) ...")
    print(f"  File name: {args.fname_sp_guesses}")

    surf.sp_guesses = []
    with open(args.fname_sp_guesses, "r", encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) >= 2:
                idx = len(surf.sp_guesses) + 1
                surf.sp_guesses.append({"id": int(idx), "x": float(fields[1])})
            else:
                raise ValueError("Stationary point guess file must contain at least 2 column.")

    print(f"  Number of loaded stationary point guesses: {len(surf.sp_guesses)}")
    print("")
    print("#       ID          X")
    print("# -------- ----------")
    for pt in surf.sp_guesses:
        print(f"{pt['id']:10d} {pt['x']:10.3f}")

    if args.showgpts or args.savegpts is not None:
        surf.plot_rbf_with_sp_guesses(show=args.showgpts, save=args.savegpts, figsize=args.figsize, dpi=args.dpi)

# ------------------------------------------------------------------------------

def sp_guess(args, surf):
    print("")
    print("# Detect stationary points ...")
    surf.detect_sp_regions_connected(eps=args.sng_cutoff, representative="minabs", min_size=args.sng_min_size)

    print(f"  Number of detected stationary points: {len(surf.sp_guesses)}")
    print("")
    print("# Found stationary points ...")
    print(f"  Saved as: {args.fname_sp_guesses}")

    with open(args.fname_sp_guesses, "w", encoding="utf-8") as fout_spg:
        fout_spg.write("# Automatically detected stationary-point guesses\n")
        fout_spg.write("#       ID          X     ENE(X)     SNG(X)       Size\n")
        fout_spg.write("# -------- ---------- ---------- ---------- ----------\n")
        for pt in surf.sp_guesses:
            fout_spg.write(
                f"{pt['id']:10d} "
                f"{pt['x']:10.3f} "
                f"{pt['ene']:10.3f} "
                f"{pt['sng']:10.3f} "
                f"{pt['size']:10d}\n"
            )

    print("")
    print("#       ID          X     ENE(X)     SNG(X)       Size")
    print("# -------- ---------- ---------- ---------- ----------")
    for pt in surf.sp_guesses:
        print(f"{pt['id']:10d} {pt['x']:10.3f} {pt['ene']:10.3f} {pt['sng']:10.3f} {pt['size']:10d}")

    if args.showgpts or args.savegpts is not None:
        surf.plot_rbf_with_sp_guesses(show=args.showgpts, save=args.savegpts, figsize=args.figsize, dpi=args.dpi)

    if args.showsng == True or args.savesng is not None:
        surf.plot_u_sng(show=args.showsng,save=args.savesng,figsize=args.figsize,dpi=args.dpi,zmax=args.sng_cutoff)

# ------------------------------------------------------------------------------

def load_sp_optimized(args, surf):
    print("")
    print("# Load optimized stationary points ...")
    print(f"  File name: {args.fname_sp_optimized}")

    surf.sp_optimized = []
    with open(args.fname_sp_optimized, "r", encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            line = line.split("#", 1)[0].strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) >= 7:
                surf.sp_optimized.append({
                    "id": int(fields[0]),
                    "x": float(fields[1]),
                    "ene": float(fields[2]),
                    "sng": float(fields[3]),
                    "r1": float(fields[4]),
                    "type": str(fields[5]),
                    "sl1": float(fields[6]),
                })
            else:
                raise ValueError("Stationary point file must contain at least 7 columns.")

    print(f"  Number of loaded optimized stationary points: {len(surf.sp_optimized)}")
    print("")
    print("#       ID          X     ENE(X)     SNG(X)         R1 Type        sL1")
    print("# -------- ---------- ---------- ---------- ---------- ---- ----------")
    for opt in surf.sp_optimized:
        print(
            f"{opt['id']:10d} {opt['x']:10.3f} {opt['ene']:10.3f} "
            f"{opt['sng']:10.3f} {opt['r1']:10.3f} {opt['type']:>4} {opt['sl1']:10.3f}"
        )

    if args.showopts or args.saveopts is not None:
        surf.plot_rbf_with_sp_optimized(show=args.showopts, save=args.saveopts, figsize=args.figsize, dpi=args.dpi)

# ------------------------------------------------------------------------------

def opt_sps(args, surf):
    print("")
    print("# Optimize stationary points ...")

    for gidx, gpt in enumerate(surf.sp_guesses,start=1):
        print("")
        print("# ==============================================================================")
        print(f">>> Processing initial guess #{gidx}: x = {gpt['x']}")
        rsize = None
        if "size" in gpt:
            rsize = gpt["size"]
        surf.find_sp(gpt["x"],trust_r0=args.trust_r0,rsize=rsize)

    print("# ==============================================================================")
    print("")
    print("# Cleaning stationary points ...")

    # sort by x
    surf.sp_optimized = sorted(surf.sp_optimized, key=lambda p: p["x"])

    print("")
    print("#       ID          X     ENE(X)     SNG(X)         R1 Type        sL1")
    print("# -------- ---------- ---------- ---------- ---------- ---- ----------")
    for opt in surf.sp_optimized:
        print(
            f"{opt['id']:10d} {opt['x']:10.3f} {opt['ene']:10.3f} "
            f"{opt['sng']:10.3f} {opt['r1']:10.3f} {opt['type']:>4} {opt['sl1']:10.3f}"
        )

    surf.remove_sng_outliers(max_sng=args.max_sng)
    surf.remove_stationary_point_outliers(minslam=args.minslam,maxslam=args.maxslam)
    surf.remove_near_duplicate_sp_optimized(min_distance_u=args.min_distance_u)

    print("")
    print("# Optimized stationary points ...")
    print(f"  Saved as:                        {args.fname_sp_optimized}")
    print(f"  Number of all stationary points: {len(surf.sp_optimized)}")
    print(f"  Number of minima:                {sum(1 for item in surf.sp_optimized if item['type'] == 'S')}")
    print(f"  Number of transition states:     {sum(1 for item in surf.sp_optimized if item['type'] == 'T')}")

    with open(args.fname_sp_optimized, "w", encoding="utf-8") as fout_spo:
        fout_spo.write("# Optimized stationary points\n")
        fout_spo.write("#       ID          X     ENE(X)     SNG(X)         R1 Type        sL1\n")
        fout_spo.write("# -------- ---------- ---------- ---------- ---------- ---- ----------\n")
        for opt in surf.sp_optimized:
            fout_spo.write(
                f"{opt['id']:10d} "
                f"{opt['x']:10.3f} "
                f"{opt['ene']:10.3f} "
                f"{opt['sng']:10.3f} "
                f"{opt['r1']:10.3f} "
                f"{opt['type']:>4s} "
                f"{opt['sl1']:10.3f}\n"
            )

    print("")
    print("#       ID          X     ENE(X)     SNG(X)         R1       Type")
    print("# -------- ---------- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        print(
            f"{opt['id']:10d} {opt['x']:10.3f} {opt['ene']:10.3f} "
            f"{opt['sng']:10.3f} {opt['r1']:10.3f} {opt['type']:>10}"
        )

    if args.showopts or args.saveopts is not None:
        surf.plot_rbf_with_sp_optimized(show=args.showopts, save=args.saveopts, figsize=args.figsize, dpi=args.dpi)

# ------------------------------------------------------------------------------

def find_basins(args, surf):
    print("")
    print("# Find minimum basins ...")

    surf.calculate_basins()
    surf.calculate_basins_ene_state()

    print(f"  Number of minimum basins:        {len(surf.basins)}")

    with open(args.fname_sp_basins, "w") as fout_spo:
        fout_spo.write("# Basin stationary points\n")
        fout_spo.write("#       ID          X     ENE(X)     SNG(X)         R1 Type\n")
        fout_spo.write("# -------- ---------- ---------- ---------- ---------- ----\n")
        for opt in surf.sp_optimized:
            fout_spo.write(
                f"{opt['id']:10d} "
                f"{opt['x']:10.3f} "
                f"{opt['ene']:10.3f} "
                f"{opt['sng']:10.3f} "
                f"{opt['r1']:10.3f} "
                f"{opt['type']:>4s}\n"
            )

    print("")
    print("#       ID Type          X       Xmin       Xmax   A=ENE(X)      A0(X)   A(state)  A0(state)")
    print("# -------- ---- ---------- ---------- ---------- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        if opt['type'] == "S":
            print(
                f"{opt['id']:10d} {opt['type']:>4} {opt['x']:10.3f} "
                f"{opt['xmin']:10.3f} {opt['xmax']:10.3f} "
                f"{opt['ene']:10.3f} {opt['ene0']:10.3f} {opt.get('ene_state', np.nan):10.3f} {opt.get('ene_state0', np.nan):10.3f}"
            )
    print("")

    with open(args.fname_basins, "w") as fout_spo:
        fout_spo.write("# Basins\n")
        fout_spo.write("#       ID Type          X       Xmin       Xmax   A=ENE(X)      A0(X)   A(state)  A0(state)\n")
        fout_spo.write("# -------- ---- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n")

        for opt in surf.sp_optimized:
            if opt['type'] == "S":
                fout_spo.write(
                    f"{opt['id']:10d} "
                    f"{opt['type']:>4s} "
                    f"{opt['x']:10.3f} "
                    f"{opt['xmin']:10.3f} {opt['xmax']:10.3f} "
                    f"{opt['ene']:10.3f} "
                    f"{opt['ene0']:10.3f} "
                    f"{opt['ene_state']:10.3f} "
                    f"{opt['ene_state0']:10.3f}\n"
                )

    if args.showbasins or args.savebasins is not None:
        surf.plot_rbf_with_basins(show=args.showbasins, save=args.savebasins, figsize=args.figsize, dpi=args.dpi)

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    print("")
    print("# ==============================================================================")
    print("#                        *** Analyze 1D Energy Surface ***                     #")
    print("#          The analyse-1D-surface utility is part of the PMFLib toolkit.       #")
    print("#==============================================================================#")
    print("# PMFLib - Potential of Mean Force Toolkit                                     #")
    print("# -----------------------------------------------------------------------------#")
    print("# Authors: (c) 2026 Petr Kulhanek (NCBR)                                       #")
    print("#                                                                              #")
    print("# NCBR:    National Centre for Biomolecular Research, Masaryk University, CZ   #")
    print("#                                                                              #")
    print("# PMFLib is licensed under Lesser GPL v2.1 and above.                          #")
    print("#==============================================================================#")

    print("")
    args = parse_args()

    options = vars(args)
    print("# All arguments ...")
    print(options)

    print("")
    print("# Initializing 1D energy surface ...")

    surf = EnergySurface1D(
        cv1min=args.cv1min,
        cv1max=args.cv1max,
        cv1_nrbfs=args.cv1nrbfs,
        cv1_nbins=args.cv1nbins,
        cv1_label=args.cv1label,
        zmax=args.zmax,
        ene_label=args.enelabel,
        temp=args.temp,
        thrfac=args.thrfac,
        random_seed=args.random_seed,
    )

    actions = [item.strip() for item in args.actions.split(",") if item.strip()]
    if "all" in actions:
        actions = ["loadfes", "optrbfs", "calcsurfs", "guess", "optimize", "basins"]

    print(f"  Actions: {', '.join(actions)}")

    for action in actions:
        if action == "loadfes":
            load_fes(args, surf)
        elif action == "optrbfs":
            opt_rbfs(args, surf)
        elif action == "calcsurfs":
            cal_surfs(args, surf)
        elif action == "loadguess":
            load_sp_guesses(args, surf)
        elif action == "guess":
            sp_guess(args, surf)
        elif action == "loadopts":
            load_sp_optimized(args, surf)
        elif action == "optimize":
            opt_sps(args, surf)
        elif action == "basins":
            find_basins(args, surf)
        else:
            raise ValueError(f"Unsupported action: {action:s}")

    print("")

# ==============================================================================

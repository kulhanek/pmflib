#!/usr/bin/env python3
"""
Simplified String Method on 2D Energy Surface

This is a Python rewrite of stm_path01.m and C++ implementation of STM in PMFLib.

The FES is represented by the same basic Gaussian RBF idea used in analyse-2D-surface.py: coordinates are scaled to
[0, 1], a tensor-product Gaussian RBF basis is fitted by SVD, and energies plus analytical gradients are evaluated from the fitted surface.

In this implementation, the metric tensor (MTC) is considered as an unit matrix.
"""
# =============================================================================

from __future__ import annotations

import argparse
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Iterable
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from scipy.spatial import cKDTree

# ------------------------------------------------------------------------------

RAD2DEG = 180.0 / np.pi

# =============================================================================
# Cubic splines for path parametrisation
# =============================================================================

class CVSplineBase:
    """Small Python analogue of the PMFLib CV spline interface."""

    def __init__(self):
        self.clear()

# ------------------------------------------------------------------------------

    def clear(self):
        self.n = -1
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)
        self.sa = np.zeros(0, dtype=float)
        self.sb = np.zeros(0, dtype=float)
        self.sc = np.zeros(0, dtype=float)
        self.sd = np.zeros(0, dtype=float)

# ------------------------------------------------------------------------------

    def allocate(self, numofknots: int):
        self.clear()
        self.n = int(numofknots) - 1
        if self.n <= 0:
            self.n = 0
            return
        size = self.n + 1
        self.x = np.zeros(size, dtype=float)
        self.y = np.zeros(size, dtype=float)
        self.sa = np.zeros(size, dtype=float)
        self.sb = np.zeros(size, dtype=float)
        self.sc = np.zeros(size, dtype=float)
        self.sd = np.zeros(size, dtype=float)

# ------------------------------------------------------------------------------

    def set_point(self, knotid: int, alpha: float, cv: float):
        if knotid < 0 or knotid > self.n:
            raise RuntimeError("knotid is out-of-range")
        self.x[knotid] = float(alpha)
        self.y[knotid] = float(cv)

# ------------------------------------------------------------------------------

    def update_points(self, alphas, values):
        alphas = np.asarray(alphas, dtype=float)
        values = np.asarray(values, dtype=float)
        if alphas.ndim != 1 or values.ndim != 1:
            raise ValueError("alphas and values must be one-dimensional arrays")
        if len(alphas) != len(values):
            raise ValueError("alphas and values must have the same length")
        if len(alphas) < 2:
            raise RuntimeError("at least two knots are required")
        if np.any(np.diff(alphas) <= 0.0):
            raise RuntimeError("spline knots must have strictly increasing alpha values")
        if self.n != len(alphas) - 1:
            self.allocate(len(alphas))
        self.x[:] = alphas
        self.y[:] = values
        self.build_spline()

# ------------------------------------------------------------------------------

    def _interval_index(self, alpha: float) -> tuple[int, float]:
        if self.n <= 0:
            raise RuntimeError("not enough of knots")
        a = float(np.clip(alpha, self.x[0], self.x[self.n]))
        i = int(np.searchsorted(self.x, a, side="right") - 1)
        i = min(max(i, 0), self.n - 1)
        return i, a - self.x[i]

# ------------------------------------------------------------------------------

    def get_cv(self, alpha: float) -> float:
        i, dx = self._interval_index(alpha)
        return float(self.sd[i] + self.sc[i]*dx + self.sb[i]*dx*dx + self.sa[i]*dx*dx*dx)

# ------------------------------------------------------------------------------

    def get_cv_first_der(self, alpha: float) -> float:
        i, dx = self._interval_index(alpha)
        return float(self.sc[i] + 2.0*self.sb[i]*dx + 3.0*self.sa[i]*dx*dx)

# ------------------------------------------------------------------------------

    def __call__(self, alpha, der: int = 0):
        if der not in (0, 1):
            raise ValueError("only value (der=0) and first derivative (der=1) are supported")
        arr = np.asarray(alpha, dtype=float)
        if arr.ndim == 0:
            return self.get_cv_first_der(float(arr)) if der == 1 else self.get_cv(float(arr))
        fn = self.get_cv_first_der if der == 1 else self.get_cv
        return np.array([fn(a) for a in arr], dtype=float)

# =============================================================================

class CVSplineInterpolatingCubic(CVSplineBase):
    """Natural interpolating cubic spline translated from CCVSplineInterpolatingCubic."""

    def build_spline(self):
        self.sa[:] = 0.0
        self.sb[:] = 0.0
        self.sc[:] = 0.0
        self.sd[:] = 0.0

        if self.n <= 0:
            raise RuntimeError("not enough of knots")

        if self.n == 1:
            dx = self.x[1] - self.x[0]
            if dx == 0.0:
                raise RuntimeError("zero spline interval")
            self.sd[0] = self.y[0]
            self.sc[0] = (self.y[1] - self.y[0]) / dx
            return

        h = np.zeros(self.n + 1, dtype=float)
        p = np.zeros(self.n + 1, dtype=float)
        q = np.zeros(self.n + 1, dtype=float)
        b = np.zeros(self.n + 1, dtype=float)

        h[0] = self.x[1] - self.x[0]
        for i in range(1, self.n):
            h[i] = self.x[i + 1] - self.x[i]
            p[i] = 2.0 * (self.x[i + 1] - self.x[i - 1])
            q[i] = 3.0*(self.y[i + 1] - self.y[i])/h[i] - 3.0*(self.y[i] - self.y[i - 1])/h[i - 1]

        for i in range(2, self.n):
            p[i] = p[i] - h[i - 1]*h[i - 1]/p[i - 1]
            q[i] = q[i] - q[i - 1]*h[i - 1]/p[i - 1]

        b[self.n - 1] = q[self.n - 1]/p[self.n - 1]
        for i in range(2, self.n):
            j = self.n - i
            b[j] = (q[j] - h[j]*b[j + 1]) / p[j]

        self.sa[0] = b[1] / (3.0*h[0])
        self.sb[0] = 0.0
        self.sc[0] = (self.y[1] - self.y[0])/h[0] - b[1]*h[0]/3.0
        self.sd[0] = self.y[0]

        for i in range(1, self.n):
            self.sa[i] = (b[i + 1] - b[i]) / (3.0*h[i])
            self.sb[i] = b[i]
            self.sc[i] = (b[i] + b[i - 1])*h[i - 1] + self.sc[i - 1]
            self.sd[i] = self.y[i]

# =============================================================================

class CVSplineSmoothingCubic(CVSplineBase):
    """Smoothing cubic spline translated from CCVSplineSmoothingCubic."""

    def __init__(self, lam: float = 0.999, sigma: float = 0.01):
        self.lambda_ = float(lam)
        self.all_sigma = float(sigma)
        super().__init__()

# ------------------------------------------------------------------------------

    def clear(self):
        super().clear()
        self.sigma = np.zeros(0, dtype=float)

# ------------------------------------------------------------------------------

    def allocate(self, numofknots: int):
        super().allocate(numofknots)
        if self.n > 0:
            self.sigma = np.full(self.n + 1, self.all_sigma, dtype=float)
        else:
            self.sigma = np.zeros(0, dtype=float)

# ------------------------------------------------------------------------------

    def set_point(self, knotid: int, alpha: float, cv: float):
        super().set_point(knotid, alpha, cv)
        self.sigma[knotid] = self.all_sigma

# ------------------------------------------------------------------------------

    def set_lambda(self, lam: float):
        lam = float(lam)
        if lam <= 0.0 or lam > 1.0:
            raise RuntimeError("lambda out-of-range (0.0;1.0>")
        self.lambda_ = lam

# ------------------------------------------------------------------------------

    def set_sigma(self, knotid: int, sig: float):
        if knotid < 0 or knotid > self.n:
            raise RuntimeError("knotid is out-of-range")
        self.sigma[knotid] = float(sig)

# ------------------------------------------------------------------------------

    def update_points(self, alphas, values):
        alphas = np.asarray(alphas, dtype=float)
        values = np.asarray(values, dtype=float)
        if self.n != len(alphas) - 1:
            self.allocate(len(alphas))
        else:
            self.sigma[:] = self.all_sigma
        super().update_points(alphas, values)

# ------------------------------------------------------------------------------

    def build_spline(self):
        self.sa[:] = 0.0
        self.sb[:] = 0.0
        self.sc[:] = 0.0
        self.sd[:] = 0.0

        if self.n <= 0:
            raise RuntimeError("not enough of knots")

        if self.n == 1:
            dx = self.x[1] - self.x[0]
            if dx == 0.0:
                raise RuntimeError("zero spline interval")
            self.sd[0] = self.y[0]
            self.sc[0] = (self.y[1] - self.y[0]) / dx
            return

        h = np.zeros(self.n + 1, dtype=float)
        r = np.zeros(self.n + 2, dtype=float)
        f = np.zeros(self.n + 2, dtype=float)
        p = np.zeros(self.n + 1, dtype=float)
        q = np.zeros(self.n + 1, dtype=float)
        u = np.zeros(self.n + 1, dtype=float)
        v = np.zeros(self.n + 1, dtype=float)
        w = np.zeros(self.n + 1, dtype=float)

        mu = 2.0 * (1.0 - self.lambda_) / (3.0 * self.lambda_)

        h[0] = self.x[1] - self.x[0]
        r[0] = 3.0 / h[0]
        for i in range(1, self.n):
            h[i] = self.x[i + 1] - self.x[i]
            r[i] = 3.0 / h[i]
            f[i] = -(r[i - 1] + r[i])
            p[i] = 2.0 * (self.x[i + 1] - self.x[i - 1])
            q[i] = 3.0*(self.y[i + 1] - self.y[i])/h[i] - 3.0*(self.y[i] - self.y[i - 1])/h[i - 1]

        v[0] = h[0]
        for i in range(1, self.n):
            u[i] = (r[i - 1]*r[i - 1]*self.sigma[i - 1]
                  + f[i]*f[i]*self.sigma[i]
                  + r[i]*r[i]*self.sigma[i + 1])
            u[i] = mu*u[i] + p[i]
            v[i] = f[i]*r[i]*self.sigma[i] + r[i]*f[i + 1]*self.sigma[i + 1]
            v[i] = mu*v[i] + h[i]
            w[i] = mu*r[i]*r[i + 1]*self.sigma[i + 1]

        self._quincunx(u, v, w, q)

        self.sd[0] = self.y[0] - mu*r[0]*q[1]*self.sigma[0]
        self.sd[1] = self.y[1] - mu*(f[1]*q[1] + r[1]*q[2])*self.sigma[0]
        self.sa[0] = q[1] / (3.0*h[0])
        self.sb[0] = 0.0
        self.sc[0] = (self.sd[1] - self.sd[0])/h[0] - q[1]*h[0]/3.0
        r[0] = 0.0

        for j in range(1, self.n):
            self.sa[j] = (q[j + 1] - q[j]) / (3.0*h[j])
            self.sb[j] = q[j]
            self.sc[j] = (q[j] + q[j - 1])*h[j - 1] + self.sc[j - 1]
            self.sd[j] = r[j - 1]*q[j - 1] + f[j]*q[j] + r[j]*q[j + 1]
            self.sd[j] = self.y[j] - mu*self.sd[j]*self.sigma[j]

# ------------------------------------------------------------------------------

    def _quincunx(self, u, v, w, q):
        u[0] = 0.0
        v[1] = v[1] / u[1]
        w[1] = w[1] / u[1]

        for j in range(2, self.n):
            u[j] = u[j] - u[j - 2]*w[j - 2]*w[j - 2] - u[j - 1]*v[j - 1]*v[j - 1]
            v[j] = (v[j] - u[j - 1]*v[j - 1]*w[j - 1]) / u[j]
            w[j] = w[j] / u[j]

        q[1] = q[1] - v[0]*q[0]
        for j in range(2, self.n):
            q[j] = q[j] - v[j - 1]*q[j - 1] - w[j - 2]*q[j - 2]

        for j in range(1, self.n):
            q[j] = q[j] / u[j]

        q[self.n] = 0.0
        for j in range(self.n - 2, 0, -1):
            q[j] = q[j] - v[j]*q[j + 1] - w[j]*q[j + 2]

# =============================================================================
# RBF surface model, adapted from analyse-2D-surface.py
# =============================================================================

# =============================================================================
# Axis
# =============================================================================

class Axis:
    """One collective variable axis, represented internally on the scaled interval [0, 1]."""

    def __init__(self, cvmin, cvmax, nrbfs, npts, label, name, stype, maxmove):
        
        if cvmax <= cvmin:
            raise ValueError("cvmax must be larger than cvmin")
        if nrbfs <= 0:
            raise ValueError("nrbfs must be positive")
        if npts <= 0:
            raise ValueError("npts must be positive")

        self.cvmin      = float(cvmin)
        self.cvmax      = float(cvmax)
        self.nrbfs      = int(nrbfs)
        self.npts       = int(npts)
        self.label      = label
        self.name       = name
        self.type       = stype
        self.maxmove    = maxmove  


        self.range = self.cvmax - self.cvmin
        self.width = 1.0 / self.nrbfs
        self.width_scale = 1.0

        if self.maxmove is None:
            self.maxmove = self.range / 10.0

        self.smaxmove = self.maxmove / self.range

        self.centers = np.arange(self.nrbfs + 1, dtype=float) / self.nrbfs

# ------------------------------------------------------------------------------

    def scale(self, x):
        u = (np.asarray(x, dtype=float) - self.cvmin) / self.range
        return u

# ------------------------------------------------------------------------------

    def unscale(self, u):
        u = np.asarray(u, dtype=float)
        return self.cvmin + u * self.range

# ------------------------------------------------------------------------------

    def delta(self, u):
        du = float(u) - self.centers
        return du

# =============================================================================
# EnergySurface2D
# =============================================================================

class EnergySurface2D:
    """2D Gaussian RBF energy surface with analytical gradients."""

    def __init__(self,x_axis,y_axis,ene_label,zmax,contour_spacing) -> None:

        self.x_axis = x_axis
        self.y_axis = y_axis
        
        self.ene_label = ene_label
        self.zmax = float(zmax)
        self.contour_spacing = float(contour_spacing)

        self.x_data: np.ndarray | None = None
        self.y_data: np.ndarray | None = None
        self.e_data: np.ndarray | None = None
        self.amplitudes: np.ndarray | None = None

        self.X: np.ndarray | None = None
        self.Y: np.ndarray | None = None
        self.Z: np.ndarray | None = None

# ------------------------------------------------------------------------------

    def load(self, filename: str | Path) -> None:
        """Load text FES data. First three columns are CV1, CV2, energy."""
        rows = []
        with open(filename, "r", encoding="utf-8") as fin:
            for line in fin:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 3:
                    continue
                rows.append([float(v) for v in fields])

        if not rows:
            raise ValueError(f"No valid FES points were loaded from {filename!s}")

        data = np.asarray(rows, dtype=float)

        self.x_data = data[:, 0]
        self.y_data = data[:, 1]
        self.e_data = data[:, 2]

        zmin = np.nanmin(self.e_data)
        if zmin < 0:
            self.e_data = self.e_data - zmin

        if self.zmax == None:
            self.zmax = np.nanmax(self.e_data)

# ------------------------------------------------------------------------------

    def to_scaled(self, x: Iterable[float]) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if x.shape != (2,):
            raise ValueError("x must be a two-element vector")
        return np.array([self.x_axis.scale(x[0]), self.y_axis.scale(x[1])], dtype=float)

# ------------------------------------------------------------------------------

    def from_scaled(self, u: Iterable[float]) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        if u.shape != (2,):
            raise ValueError("u must be a two-element vector")
        return np.array([self.x_axis.unscale(u[0]), self.y_axis.unscale(u[1])], dtype=float)

# ------------------------------------------------------------------------------

    def _basis_1d(self, u: float, axis: Axis):
        du = axis.delta(u)
        s = axis.width * axis.width_scale
        g = np.exp(-0.5 * (du / s) ** 2)
        dg = g * (-du / s**2)
        d2g = g * ((du**2 / s**4) - (1.0 / s**2))
        return g, dg, d2g

# ------------------------------------------------------------------------------

    def _basis_2d_uv(self, u: Iterable[float], derivatives: bool = False):
        u = np.asarray(u, dtype=float)
        gx, dgx, d2gx = self._basis_1d(u[0], self.x_axis)
        gy, dgy, d2gy = self._basis_1d(u[1], self.y_axis)

        phi = np.outer(gx, gy)
        if not derivatives:
            return phi.ravel()

        phi_u = np.outer(dgx, gy)
        phi_v = np.outer(gx, dgy)
        phi_uu = np.outer(d2gx, gy)
        phi_vv = np.outer(gx, d2gy)
        phi_uv = np.outer(dgx, dgy)
        return phi.ravel(), phi_u.ravel(), phi_v.ravel(), phi_uu.ravel(), phi_vv.ravel(), phi_uv.ravel()

# ------------------------------------------------------------------------------

    def _build_design_matrix(self) -> np.ndarray:
        if self.x_data is None or self.y_data is None:
            raise RuntimeError("No FES data loaded")
        nbasis = len(self.x_axis.centers) * len(self.y_axis.centers)
        B = np.empty((len(self.x_data), nbasis), dtype=float)
        for k, (x, y) in enumerate(zip(self.x_data, self.y_data)):
            B[k, :] = self._basis_2d_uv(self.to_scaled([x, y]), derivatives=False)
        return B
# ------------------------------------------------------------------------------

    def fit(self, sx=1.0, sy=1.0, rcond=1.0e-12) -> None:
        if self.e_data is None:
            raise RuntimeError("No FES data loaded")

        self.x_axis.width_scale = float(sx)
        self.y_axis.width_scale = float(sy)

        B = self._build_design_matrix()
        U, s, Vt = np.linalg.svd(B, full_matrices=False)
        cutoff = float(rcond) * np.max(s)
        sinv = np.zeros_like(s)
        keep = s > cutoff
        sinv[keep] = 1.0 / s[keep]
        coeff = Vt.T @ (sinv * (U.T @ self.e_data))
        self.amplitudes = coeff.reshape((len(self.x_axis.centers), len(self.y_axis.centers)))

        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

# ------------------------------------------------------------------------------

    def eval_uv(self, u: Iterable[float]):
        if self.amplitudes is None:
            raise RuntimeError("The RBF surface has not been fitted")
        A = self.amplitudes.ravel()
        phi, phi_u, phi_v, phi_uu, phi_vv, phi_uv = self._basis_2d_uv(u, derivatives=True)
        value = float(np.dot(A, phi))
        gradient = np.array([np.dot(A, phi_u), np.dot(A, phi_v)], dtype=float)
        hessian = np.array(
            [[np.dot(A, phi_uu), np.dot(A, phi_uv)], [np.dot(A, phi_uv), np.dot(A, phi_vv)]],
            dtype=float,
        )
        return value, gradient, hessian

# ------------------------------------------------------------------------------

    def eval(self, x: Iterable[float]):
        value = self.eval_uv(self.to_scaled(x))[0]
        return value

# ------------------------------------------------------------------------------

    def calc_grid(self) -> None:

        x_grid = np.linspace(self.x_axis.cvmin, self.x_axis.cvmax, self.x_axis.npts+1)
        y_grid = np.linspace(self.y_axis.cvmin, self.y_axis.cvmax, self.y_axis.npts+1)
        self.X, self.Y = np.meshgrid(x_grid, y_grid, indexing="xy")
        self.Z = np.empty_like(self.X, dtype=float)
        for i in range(self.Z.shape[0]):
            for j in range(self.Z.shape[1]):
                self.Z[i, j] = self.eval([self.X[i, j], self.Y[i, j]])

        # detect unsampled points
        self.unsampled_mask, nearest_dist = self.detect_unsampled_grid_points()

        # set them to zmax
        self.Z[self.unsampled_mask] = self.zmax

# ------------------------------------------------------------------------------

    def detect_unsampled_grid_points(self, max_distance=None):
        """
        Detect grid points that are too far from any sampled scattered data point.

        Parameters
        ----------
        max_distance : float or None
            Maximum allowed distance to the nearest sampled point.
            If None, it is estimated from the regular grid spacing.

        Returns
        -------
        unsampled_mask : ndarray, shape like self.X
            True where the grid point is considered unsampled.

        nearest_dist : ndarray, shape like self.X
            Distance to the nearest sampled scattered point.
        """

        if self.x_data is None or self.y_data is None:
            raise RuntimeError("No scattered data loaded. Call load() first.")

        if self.X is None or self.Y is None:
            raise RuntimeError("Regular grid is not available. Call calc_z_and_sp() first.")

        # ------------------------------------------------------------
        # Work either in scaled coordinates
        # ------------------------------------------------------------

        data_points = np.array(
            [self.to_scaled([x, y]) for x, y in zip(self.x_data, self.y_data)]
        )

        grid_points = np.column_stack([
            self.x_axis.scale(self.X.ravel()),
            self.y_axis.scale(self.Y.ravel()),
        ])

        dx = 1.0 / (self.x_axis.npts - 1)
        dy = 1.0 / (self.y_axis.npts - 1)

        # ------------------------------------------------------------
        # Default threshold:
        # a grid point is sampled if there is a data point roughly within
        # one grid-cell diagonal.
        # ------------------------------------------------------------

        if max_distance is None:
            max_distance = 0.75 * np.sqrt(dx**2 + dy**2)

        tree = cKDTree(data_points)
        nearest_dist, nearest_idx = tree.query(grid_points, k=1)

        nearest_dist = nearest_dist.reshape(self.X.shape)
        unsampled_mask = nearest_dist > max_distance

        return unsampled_mask, nearest_dist
    
# ------------------------------------------------------------------------------

    def make_colormap(self):
        colors = ["#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F"]
        bounds = np.arange(0.0, self.zmax + self.contour_spacing, self.contour_spacing)
        cmap = LinearSegmentedColormap.from_list("gnuplot_like_fes", colors, len(bounds))
        norm = BoundaryNorm(bounds, len(bounds), clip=True)
        return cmap, norm
    
# ------------------------------------------------------------------------------

    def plot_fes_with_path(self, title: str = None, 
                           path_x: np.ndarray = None, path_y: np.ndarray = None, 
                           filename: str=None, show: bool = None, dpi: int = 300) -> None:

        cmap, norm = self.make_colormap()
        cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing
        levels = np.arange(0.0, cmax + self.contour_spacing, self.contour_spacing)

        fig, ax = plt.subplots()
        im = ax.pcolormesh(self.X, self.Y, self.Z, shading="auto", cmap=cmap, norm=norm)
        fig.colorbar(im, ax=ax, label=self.ene_label)
        ax.contour(self.X, self.Y, self.Z, levels=levels, colors="black", linewidths=0.5)

        if( path_x is not None and path_y is not None):
            ax.plot(path_x, path_y, ".-", color="white", markersize=5, linewidth=1.0)
        
        ax.set_xlabel(self.x_axis.label)
        ax.set_ylabel(self.y_axis.label)
        ax.set_title(title)

        ax.set_aspect("equal", adjustable="box")
        fig.tight_layout()

        if filename is not None:
            fig.savefig(filename, dpi=dpi)

        if show:
            plt.show()
        
        plt.close(fig)

# =============================================================================
# String-method utilities
# =============================================================================

# =============================================================================
# Bead
# =============================================================================

class Bead:
    def __init__(self, stmpath: STMPath):

        self.ncvs       = stmpath.ncvs

        self.Pos        = np.zeros(self.ncvs)               #  bead position, scaled
        self.dCVdAlpha  = np.zeros(self.ncvs)
        self.P          = np.zeros((self.ncvs,self.ncvs))   # projector

        self.Grad         = np.zeros(self.ncvs)
        self.pGrad        = np.zeros(self.ncvs)             # force acting perpendicularly to the path

        self.Alpha      = None      # path position
        self.dAdAlpha   = None      # free energy derivative
        self.A          = None      # free energy, integrated
        self.Asurf      = None      # free energy from EnergySurface2D
             
        # Adam (Adaptive Moment Estimation) variants
        self.beta1t    = 1.0
        self.beta2t    = 1.0  
        self.mt        = np.zeros(self.ncvs)
        self.vt        = np.zeros(self.ncvs)

        # old pos
        self.OPos      = np.zeros(self.ncvs)

# ------------------------------------------------------------------------------

    def UpdatePositionADABelif(self,step,beta1,beta2,mingnormeps,cvs):

        self.mt[:] = beta1 * self.mt[:] + (1.0 - beta1) * self.pGrad[:]
        self.vt[:] = beta2 * self.vt[:] + (1.0 - beta2) * ((self.pGrad[:] - self.mt[:])*(self.pGrad[:] - self.mt[:]) + mingnormeps)

        self.beta1t = self.beta1t * beta1
        self.beta2t = self.beta2t * beta2

        for i in range(self.ncvs):

            mthat = self.mt[i]/(1.0-self.beta1t)
            vthat = self.vt[i]/(1.0-self.beta2t)

            dm = step * mthat/(math.sqrt(vthat)+mingnormeps)

            if (cvs[i].smaxmove <= 0) or (math.fabs(dm) < cvs[i].smaxmove):
                self.Pos[i] = self.Pos[i] - dm
            else:
                self.Pos[i] = self.Pos[i] - cvs[i].smaxmove*math.copysign(1.0,dm)

# =============================================================================
# STMPath
# =============================================================================

class STMPath:
    def __init__(self, args):

        # path parameters
        self.PathName               = args.pathname
        self.ncvs                   = 2
        self.nbeads                 = args.nbeads
        self.freeterminals          = args.freeterminals    # are path ends permanent?
        self.segment_discretization = args.segdisc          # path segment discretization

        # axes/CVs
        self.cvs = []

        x_axis = Axis(args.cv1min,args.cv1max,args.cv1nrbfs,args.cv1nbins,args.cv1label,args.cv1name,args.cv1type,args.cv1maxmove)
        self.cvs.append(x_axis)

        y_axis = Axis(args.cv2min,args.cv2max,args.cv2nrbfs,args.cv2nbins,args.cv2label,args.cv2name,args.cv2type,args.cv2maxmove)
        self.cvs.append(y_axis)

        self.surface = EnergySurface2D(x_axis,y_axis,args.enen,args.zmax,args.contour_spacing)

        print("")
        print(f"# Load FES: {args.input}")
        self.surface.load(args.input)

        print(f"")
        print(f"# Optimize RBF ...")
        print(f"  Sx:           {args.sx:10.3f}")
        print(f"  Sy:           {args.sy:10.3f}")

        self.surface.fit(sx=args.sx, sy=args.sy, rcond=args.rcond)
        print(f"  RBF fit RMSE: {self.surface.fit_rmse:10.3f}")

        print("")
        print("# Calculate Z surface ...")
        self.surface.calc_grid()

        print("")
        print("# CV splines ...")
        if args.cvspline == 0:
            print(f"  >>> Interpolating Cubic Spline")
            self.cv_splines =  [CVSplineInterpolatingCubic() for _ in range(self.ncvs)]
        else:
            print(f"  >>> Smoothing Cubic Spline")
            print(f"      Lambda: {args.spline_lambda:10.5f}")
            print(f"      Sigma:  {args.spline_sigma:10.5f}")
            self.cv_splines =  [
                CVSplineSmoothingCubic(lam=args.spline_lambda, sigma=args.spline_sigma)
                for _ in range(self.ncvs)
            ]

        print("")
        print("# Initial path ...")

        # read initial path
        input_beads = self.create_initial_beads_from_args(args)

        # plot the user initial path
        if args.plot :
            path_x, path_y = self.beads_to_xy_arrays(input_beads)
            self.surface.plot_fes_with_path(title="Initial Path - User Input",
                path_x=path_x, path_y=path_y,
                filename=f"{args.plot_prefix}_0000_a_initial_path-user.png", show=args.show, dpi=args.dpi)

        #  generate completed path
        self.beads = self.generate_beads_from_input_beads(
                        input_beads=input_beads,
                        num_beads=self.nbeads,
                    )

        # plot the user completed path
        if args.plot :
            path_x, path_y = self.beads_to_xy_arrays(self.beads)
            self.surface.plot_fes_with_path(title="Initial Path - Full Path",
                path_x=path_x, path_y=path_y,
                filename=f"{args.plot_prefix}_0000_b_initial_path-full.png", show=args.show, dpi=args.dpi)

        # setup ADABelif
        self.StepSize           = args.stepsize
        self.AdamB1             = args.beta1
        self.AdamB2             = args.beta2
        self.MinGNormEps        = args.mingnormesp

        # smoothing, reparameterization
        self.SmoothInterval     = args.smoothinterval
        self.SmoothingFac       = args.smoothingfac
        self.ReparamInterval    = args.reparaminterval

        # STM
        self.STMStep            = 0
        self.nstepmax           = args.nstepmax

        # run initial statistics
        self.init_stm_statistics(args)
        
# ------------------------------------------------------------------------------

    def stm_optimize(self,args):

        print("")
        print("# STM path optimization ...")

        self.calc_beads()
        self.integrate_path()
        self.calculate_stm_step_stat()

        if args.output is not None:
            fout=open(args.output,"w")
            self.print_path(fout=fout,beads=self.beads)
            fout.close()

        print("")
        print("# Initial path summary ...")
        self.print_path_summary_header()
        self.print_path_summary_data()

        print("")
        self.print_stm_header_f()
        self.print_stm_step_info_f()

        if args.trajectory is not None:
            ftraj = open(args.trajectory,"w")
            self.write_trajectory_header(ftraj)
            self.write_trajectory_snapshot(ftraj)

        if args.optlog is not None:
            flog=open(args.optlog,"w")
            self.print_stm_header_f(fout=flog)
            self.print_stm_step_info_f(fout=flog)            

        for step in range(1,self.nstepmax):
            # update bead positions
            self.update_all_positions()
            self.smooth_all_positions()
            self.reparametrize_all_positions()
            self.check_boundaries_of_beads(self.beads)

            # re-evaluate path
            self.calc_beads()
            self.integrate_path()
            self.calculate_stm_step_stat()

            self.print_stm_step_info_f()
            if args.optlog is not None:
                self.print_stm_step_info_f(fout=flog)

            if args.trajectory is not None:
                self.write_trajectory_snapshot(ftraj)

            if args.plot :
                path_x, path_y = self.beads_to_xy_arrays(self.beads)
                self.surface.plot_fes_with_path(title=f"Intermediate Path #{self.STMStep:04d}",
                    path_x=path_x, path_y=path_y,
                    filename=f"{args.plot_prefix}_{self.STMStep:04d}_c_path.png", show=args.show, dpi=args.dpi)

            if self.TermCrit == 5:
                break

        if self.TermCrit == 5:
            print(">>> STM path optimization was sucessfull.")
        else:
            print(">>> STM path optimization failed.")

        print("")
        print("# Final path summary ...")
        self.print_path_summary_header()
        self.print_path_summary_data()

        if args.plot :
            path_x, path_y = self.beads_to_xy_arrays(self.beads)
            self.surface.plot_fes_with_path(title="Final Path",
                path_x=path_x, path_y=path_y,
                filename=f"{args.plot_prefix}_{self.STMStep:04d}_d_final_path.png", show=args.show, dpi=args.dpi)

        if args.output is not None:
            fout=open(args.output,"w")
            self.print_path(fout=fout,beads=self.beads)
            fout.close()

        if args.optlog is not None:
            flog.close()

        if args.trajectory is not None:
            ftraj.close()

# ------------------------------------------------------------------------------

    def update_all_positions(self):

        # backup old positions
        for bead in self.beads:
            bead.OPos[:] = bead.Pos[:]

        self.OldCPathLength = self.CurrentPathLength

        self.STMStep += 1

        # update positions - for terminals
        if self.freeterminals == True:
            self.beads[0].UpdatePositionADABelif(self.StepSize,self.AdamB1,self.AdamB2,self.MinGNormEps,self.cvs)
            self.beads[-1].UpdatePositionADABelif(self.StepSize,self.AdamB1,self.AdamB2,self.MinGNormEps,self.cvs)
        
        # for the rest of the path
        for bead in self.beads[1:-1]:
            bead.UpdatePositionADABelif(self.StepSize,self.AdamB1,self.AdamB2,self.MinGNormEps,self.cvs)

# ------------------------------------------------------------------------------

    def smooth_all_positions(self):

        if (self.SmoothInterval == 0) or (self.STMStep % self.SmoothInterval != 0 ):
            return;

        for i in range(1,self.nbeads-1):
            self.beads[i].Pos = (1.0-self.SmoothingFac)*self.beads[i].Pos + 0.5*self.SmoothingFac * (self.beads[i-1].Pos + self.beads[i+1].Pos)

# ------------------------------------------------------------------------------

    def reparametrize_all_positions(self):

        if (self.ReparamInterval == 0) or (self.STMStep % self.ReparamInterval != 0 ):
            return;

        self.optimize_path(self.beads)

        # generate evenly distributed set of alphas
        for i in range(self.nbeads):
            self.beads[i].Alpha = float(i) / float(self.nbeads-1)

        self.beads[0].Alpha = 0.0
        self.beads[-1].Alpha = 1.0

        alphas = np.array([bead.Alpha for bead in self.beads], dtype=float)
        #print(alphas)

        # update positions along the path based on new alphas
        for cv in range(self.ncvs):
            # if free terminals, update their positions
            if self.freeterminals == True:
                self.beads[0].Pos[cv]  = self.cv_splines[cv](0.0)
                self.beads[-1].Pos[cv] = self.cv_splines[cv](1.0)
            
            # redistribute other beads
            for bead in self.beads[1:-1]:
                bead.Pos[cv] = self.cv_splines[cv](bead.Alpha)

# ------------------------------------------------------------------------------

    def optimize_path(self, beads):
        """
        Determine Alpha values according to bead Positions.

        Parameters
        ----------
        beads : list[Bead]
            Beads defining the path. Each bead stores bead.Pos.

        Returns
        -------
        float
            Final spline path length.
        """

        if len(beads) < 2:
            raise RuntimeError("len(beads) must be greater or equal to 2")

        # ---------------------------------------------------------------------
        # Initial path length from linear interpolation.
        # ---------------------------------------------------------------------

        total_length = 0.0

        for b in range(1, len(beads)):
            diff = beads[b].Pos - beads[b - 1].Pos
            total_length += np.linalg.norm(diff)

        if total_length == 0.0:
            raise RuntimeError("path has zero length")

        # ---------------------------------------------------------------------
        # Initial alpha values from linear interpolation.
        # ---------------------------------------------------------------------

        beads[0].Alpha = 0.0

        path_length = 0.0

        for b in range(1, len(beads) - 1):
            diff = beads[b].Pos - beads[b - 1].Pos
            segment_length = np.linalg.norm(diff)

            if segment_length == 0.0:
                raise RuntimeError("path segment has zero length")

            path_length += segment_length
            beads[b].Alpha = path_length / total_length

        beads[-1].Alpha = 1.0

        # ---------------------------------------------------------------------
        # Iteratively improve the alpha values using spline arc length.
        # ---------------------------------------------------------------------

        for iter in range(1000):

            # Build one spline for each CV coordinate:
            #
            #   CV_i = CV_i(alpha)
            #
            alphas = np.array([bead.Alpha for bead in beads], dtype=float)

            for i in range(self.ncvs):
                values = np.array([bead.Pos[i] for bead in beads], dtype=float)
                self.cv_splines[i].update_points(alphas, values)

            previous_length = total_length

            # Determine current spline path length.
            total_length = 0.0

            for b in range(1, len(beads)):
                total_length += self.get_segment_length(
                    beads[b - 1].Alpha,
                    beads[b].Alpha,
                )

            if abs(total_length - previous_length) < 1.0e-7:
                return total_length

            # Determine new alpha values.
            beads[0].Alpha = 0.0

            path_length = 0.0
            previous_alpha = beads[0].Alpha

            for b in range(1, len(beads) - 1):
                path_length += self.get_segment_length(
                    previous_alpha,
                    beads[b].Alpha,
                )

                beads[b].Alpha = path_length / total_length
                previous_alpha = beads[b].Alpha

            beads[-1].Alpha = 1.0

        return total_length

# -------------------------------------------------------------------------

    def get_segment_length(self, alpha1, alpha2):
        alphas = np.linspace(
            alpha1,
            alpha2,
            self.segment_discretization + 1,
        )

        points = np.column_stack([
            self.cv_splines[i](alphas)
            for i in range(self.ncvs)
        ])

        diffs = points[1:] - points[:-1]
        return np.sum(np.linalg.norm(diffs, axis=1))


# -------------------------------------------------------------------------

    def create_initial_beads_from_args(self, args):
        """
        Create initial STM beads from command-line arguments.

        The input coordinates are expected in physical EnergySurface2D
        coordinates, for example:

            --initial-cv1 "-3.0,-2.5,-0.5,1.0,3.2"
            --initial-cv2 "-2.0,-1.5,-1.5,0.0,2.0"

        The returned bead positions are stored in scaled EnergySurface2D
        coordinates:

            u = surface.to_scaled([cv1, cv2])

        The bead alpha values are intentionally left uninitialized.
        They are later assigned by optimize_path().
        """

        cv1_values = args.initial_cv1
        cv2_values = args.initial_cv2

        if len(cv1_values) != len(cv2_values):
            raise ValueError(
                "--initial-cv1 and --initial-cv2 must contain the same number "
                f"of values, got {len(cv1_values)} and {len(cv2_values)}"
            )

        if len(cv1_values) < 2:
            raise ValueError("At least two initial points are required")

        beads = []

        for cv1, cv2 in zip(cv1_values, cv2_values):
            pos_scaled = self.surface.to_scaled([cv1, cv2])
            bead = Bead(self)
            bead.Pos = pos_scaled
            beads.append(bead)

        return beads
    
# -------------------------------------------------------------------------

    def generate_beads_from_input_beads(
        self,
        input_beads,
        num_beads,
        check_boundaries=True
    ):
        """
        Generate the final uniformly distributed STM bead list from user input beads.

        This is the Python analogue of the C++ setup logic:
        
            2. Optimize/reparametrize the input path.
            3. Generate NumOfBeads beads from the spline path.
            4. Check/correct boundaries.
            5. Re-optimize the corrected path.
            6. Generate final uniformly spaced bead positions.

        Notes
        -----
        All coordinates are assumed to be in scaled EnergySurface2D coordinates,
        usually [0, 1] x [0, 1].

        Parameters
        ----------
        input_beads : list[Bead]
            User-defined beads, usually created from --initial-cv1 and --initial-cv2.

        num_beads : int
            Total number of STM beads to generate.

        check_boundaries : bool
            If True, bead positions are passed through boundary correction.

        Returns
        -------
        beads : list[Bead]
            Final generated bead list.
        """

        if len(input_beads) < 2:
            raise RuntimeError("At least two input beads are required")

        if num_beads < 2:
            raise RuntimeError("num_beads must be greater or equal to 2")

        # ---------------------------------------------------------------------
        # 1) Optimize the user-provided path.
        # ---------------------------------------------------------------------

        self.optimize_path(input_beads)

        # ---------------------------------------------------------------------
        # 2) Generate missing beads from the optimized input spline.
        # ---------------------------------------------------------------------

        beads = [Bead(self) for _ in range(num_beads)]

        for b, bead in enumerate(beads):
            alpha = float(b) / float(num_beads - 1)

            bead.Alpha = alpha

            for i in range(self.ncvs):
                bead.Pos[i] = self.cv_splines[i](alpha)

        # Force exact endpoint alphas.
        beads[0].Alpha = 0.0
        beads[-1].Alpha = 1.0

        # ---------------------------------------------------------------------
        # 3) Check boundaries.
        # ---------------------------------------------------------------------

        if check_boundaries:
            self.check_boundaries_of_beads(beads)

        # ---------------------------------------------------------------------
        # 4) Re-optimize path after boundary correction.
        # ---------------------------------------------------------------------

        self.optimize_path(beads)

        # ---------------------------------------------------------------------
        # 5) Final correction: regenerate bead positions at uniform alpha values.
        # ---------------------------------------------------------------------

        for b, bead in enumerate(beads):
            alpha = float(b) / float(num_beads - 1)

            bead.Alpha = alpha

            for i in range(self.ncvs):
                bead.Pos[i] = self.cv_splines[i](alpha)

        beads[0].Alpha = 0.0
        beads[-1].Alpha = 1.0

        return beads
    
# -------------------------------------------------------------------------

    def check_boundaries_of_beads(self, beads):
        """
        Correct bead.Pos positions after path generation.

        This simple version assumes all CVs are already in scaled coordinates.
        Thus, valid coordinates are clipped to [0, 1].
        """

        for bead in beads:
            bead.Pos[:] = np.clip(bead.Pos, 0.0, 1.0)
        
# -------------------------------------------------------------------------

    def beads_to_xy_arrays(self, beads):
        if len(beads) == 0:
            raise ValueError("bead list must not be empty")

        scaled_positions = np.array([bead.Pos for bead in beads], dtype=float)

        if scaled_positions.ndim != 2 or scaled_positions.shape[1] != 2:
            raise ValueError("all bead positions must have shape (2,)")

        xy = np.array(
            [self.surface.from_scaled(pos) for pos in scaled_positions],
            dtype=float,
        )

        return xy[:, 0], xy[:, 1]

# -------------------------------------------------------------------------

    def calc_beads(self):
        """
        For all beads:
            Calculate bead Asurf and MF.
            Calculate perpendicular projectors and projected mean forces
        
        Assumptions
        -----------
        - bead.Pos, bead.Grad, bead.pGrad are in scaled coordinates.
        - MTZ is the identity matrix and is omitted.
        - No additional CV range scaling is applied.
        - self.cv_splines were already built by optimize_path()
        """

        self.CurrentPathLength = self.optimize_path(self.beads)

        for bead in self.beads:            
            # ---------------------------------------------------------------------
            # Calculate MF and A
            # ---------------------------------------------------------------------

            bead.Asurf, grad, hessian = self.surface.eval_uv(bead.Pos)
            bead.Grad[:] = grad[:]

            # ---------------------------------------------------------------------
            # Calculate tangent dCV/dalpha from the path splines.
            # ---------------------------------------------------------------------

            for i in range(self.ncvs):
                # scipy.interpolate.CubicSpline:
                # spline(alpha, 1) gives the first derivative.
                bead.dCVdAlpha[i] = self.cv_splines[i](bead.Alpha, 1)


            slen2 = float(np.dot(bead.dCVdAlpha, bead.dCVdAlpha))

            if slen2 == 0.0:
                raise RuntimeError("derivative segment has zero length")

            # ---------------------------------------------------------------------
            # Projector perpendicular to the path:
            #
            #     P = I - t t^T / |t|^2
            # ---------------------------------------------------------------------

            bead.P[:, :] = np.eye(self.ncvs) - np.outer(bead.dCVdAlpha, bead.dCVdAlpha) / slen2

        # ---------------------------------------------------------------------
        # Project mean force.
        #
        # End beads move by steepest descent, without perpendicular projection.
        # Internal beads move only perpendicular to the path.
        # ---------------------------------------------------------------------

        if( self.freeterminals ):
            self.beads[0].pGrad[:] = self.beads[0].Grad[:]
            self.beads[-1].pGrad[:] = self.beads[-1].Grad[:]
        else:
            self.beads[0].pGrad[:] = 0.0
            self.beads[-1].pGrad[:] = 0.0
        
        for bead in self.beads[1:-1]:
            bead.pGrad[:] = bead.P @ bead.Grad

# -------------------------------------------------------------------------

    def integrate_path(self):
        """
        Integrate the mean-force projection along the STM path.
        Assumptions
        -----------
        - bead.Alpha is already assigned.
        - self.cv_splines are already built.
        - bead.Grad is in scaled coordinates.
        - bead.dCVdAlpha is in scaled coordinates.
        - No CV range scaling is applied.
        - The free-energy profile A(alpha) is obtained by trapezoidal integration:

            dA/dalpha = dot(dCV/dalpha, MF)

        Parameters
        ----------
        beads : list[Bead]
            STM bead list.

        Returns
        -------
        None
            The method updates each bead in place:
                bead.dAdAlpha
                bead.A
        """

        fes = 0.0

        for b, bead in enumerate(self.beads):

            # ---------------------------------------------------------------------
            # Projection of mean force along the path.
            # Here all quantities are already in scaled coordinates, so:
            #
            #     dA/dalpha = dot(dCV/dalpha, MF)
            # ---------------------------------------------------------------------

            bead.dAdAlpha = float(np.dot(bead.dCVdAlpha, bead.Grad))

            # ---------------------------------------------------------------------
            # Trapezoidal integration along alpha.
            # ---------------------------------------------------------------------

            if b > 0:
                prev = self.beads[b - 1]

                da = bead.Alpha - prev.Alpha

                if da < 0.0:
                    raise RuntimeError("bead Alpha values must be non-decreasing")

                fes += 0.5 * da * (bead.dAdAlpha + prev.dAdAlpha)

            bead.A = fes

        # -------------------------------------------------------------------------
        # Shift global minimum to zero.
        # -------------------------------------------------------------------------

        amin = min(bead.A for bead in self.beads)

        for bead in self.beads:
            bead.A -= amin

# -------------------------------------------------------------------------

    def print_path(self, fout=None, beads=None):
            """
            Print the STM path in the same format as CSTMPath::PrintPath().

            Notes
            -----
            - Only flexible terminal beads are supported.
            - All bead positions are stored internally in scaled coordinates.
            - Printed bead positions are converted back to physical CV units.
            """

            import sys

            if fout is None:
                fout = sys.stdout

            if beads is None:
                beads = self.beads

            print("[PATH]", file=fout)
            print(f"name      {self.PathName}", file=fout)
            print(f"ncvs      {self.ncvs}", file=fout)
            print(f"nbeads    {self.nbeads}", file=fout)

            print("names     ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.name:>12}", end="", file=fout)
            print(file=fout)

            print("types     ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.type:>12}", end="", file=fout)
            print(file=fout)

            print("min       ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.cvmin:12.5e}", end="", file=fout)
            print(file=fout)

            print("max       ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.cvmax:12.5e}", end="", file=fout)
            print(file=fout)

            print("maxmov    ", end="", file=fout)
            for cv in self.cvs:
                if cv.maxmove is None:
                    print(f" {'None':>12}", end="", file=fout)
                else:
                    print(f" {cv.maxmove:12.5e}", end="", file=fout)
            print(file=fout)

            for index, bead in enumerate(beads):
                # Python version supports only flexible beads.
                if self.is_bead_permanent(index):
                    print("permanent ", end="", file=fout)
                else:
                    print("flexible  ", end="", file=fout)

                for i, cv in enumerate(self.cvs):
                    scaled = bead.Pos[i]
                    unscaled = cv.unscale(scaled)
                    print(f" {float(unscaled):12.5e}", end="", file=fout)

                print(file=fout)

# -------------------------------------------------------------------------

    def is_bead_permanent(self, bead_index):
        """
        Return True if the bead should be reported as permanent.

        The simplified Python STM implementation does not store the C++
        Bead::Permanent flag explicitly.  Fixed terminal beads are therefore
        interpreted as permanent when --freeterminals is not active.
        """

        return (not self.freeterminals) and (bead_index == 0 or bead_index == self.nbeads - 1)

# -------------------------------------------------------------------------

    def print_path_summary_header(self, fout=None):
        """
        Print the path-summary header.

        This is the Python analogue of CSTMPath::PrintPathSummaryHeader().
        All bead positions and derivatives are stored internally in scaled CV
        units.  The CV position columns are printed in physical CV units, as in
        print_path().
        """

        import sys

        if fout is None:
            fout = sys.stdout

        print("# === [PATH] ===================================================================", file=fout)
        print(f"# Path name       = {self.PathName}", file=fout)
        print(f"# Number of CVs   = {self.ncvs}", file=fout)
        print(f"# Number of beads = {self.nbeads}", file=fout)

        # Header legends.
        print("#  ID   Type  MO ST  alpha  dA/dalpha       A            CID Updates", end="", file=fout)
        for i in range(self.ncvs):
            print(f"     CV{i + 1:<2d}    ", end="", file=fout)
        for i in range(self.ncvs):
            print(f"   dA/dCV{i + 1:<2d}  ", end="", file=fout)
        for i in range(self.ncvs):
            print(f" dCV{i + 1:<2d}/dalpha", end="", file=fout)
        for i in range(self.ncvs):
            print(f" -|F{i + 1:<2d}/dalpha", end="", file=fout)
        print(file=fout)

        # Delimiters.
        delimiter = "# ---- ------ -- -- ------ ------------ ------------ ------- -------"
        delimiter += " ------------" * (4 * self.ncvs)
        print(delimiter, file=fout)

        # CV metadata rows.
        print(f"{'#      names':<68}", end="", file=fout)
        for _ in range(4):
            for cv in self.cvs:
                print(f" {cv.name:>12}", end="", file=fout)
        print(file=fout)

        print(f"{'#      types':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.type:>12}", end="", file=fout)
        print(file=fout)

        print(f"{'#      min':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.cvmin:12.5e}", end="", file=fout)
        print(file=fout)

        print(f"{'#      max':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.cvmax:12.5e}", end="", file=fout)
        print(file=fout)

        print(f"{'#      maxmov':<68}", end="", file=fout)
        for cv in self.cvs:
            if cv.maxmove is not None and cv.maxmove > 0.0:
                print(f" {cv.maxmove:12.5e}", end="", file=fout)
            else:
                print(f" {'--':>12}", end="", file=fout)
        print(file=fout)

        print(delimiter, file=fout)

        print("#    1      2  3  4      5            6            7       8       9", end="", file=fout)
        column_id = 10
        for _ in range(4):
            for _ in range(self.ncvs):
                print(f"{column_id:13d}", end="", file=fout)
                column_id += 1
        print(file=fout)

        print(delimiter, file=fout)

# -------------------------------------------------------------------------

    def print_path_summary_data(self, fout=None):
        """
        Print path-summary data rows.

        Python-field mapping relative to the C++ STM implementation:
            Beads[b]->MF  -> bead.Grad
            Beads[b]->pMF -> bead.pGrad

        The simplified standalone implementation does not have client modes,
        mode statuses, or client IDs.  These columns are therefore printed as
        '--' while preserving the original table layout.
        """

        if fout is None:
            fout = sys.stdout

        for bead_index, bead in enumerate(self.beads):
            bead_id = bead_index + 1
            bead_type = "P" if self.is_bead_permanent(bead_index) else "F"
            mode = "--"
            status = "--"
            client_id = "--"

            alpha = 0.0 if bead.Alpha is None else float(bead.Alpha)
            d_ad_alpha = 0.0 if bead.dAdAlpha is None else float(bead.dAdAlpha)
            free_energy = 0.0 if bead.A is None else float(bead.A)

            print(
                f"  {bead_id:4d} {bead_type:>6} {mode:>2} {status:>2} "
                f"{alpha:6.4f} "
                f"{d_ad_alpha:12.5e} "
                f"{free_energy:12.5e} "
                f"{client_id:>7}"
                f"{self.STMStep:8d}",
                end="",
                file=fout,
            )

            # CV coordinates in physical units.
            for i, cv in enumerate(self.cvs):
                unscaled = cv.unscale(bead.Pos[i])
                print(f" {float(unscaled):12.5e}", end="", file=fout)

            # dA/dCV in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.Grad[i]):12.5e}", end="", file=fout)

            # dCV/dalpha in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.dCVdAlpha[i]):12.5e}", end="", file=fout)

            # Projected force, i.e. perpendicular component in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.pGrad[i]):12.5e}", end="", file=fout)

            print(file=fout)

# -------------------------------------------------------------------------

    def write_trajectory_header(self, fout=None):
        
        if fout is None:
            fout = sys.stdout

        print(f"# STMTRAJ {self.ncvs} {self.nbeads}", file=fout)
        self.print_path_summary_header(fout)

# -------------------------------------------------------------------------

    def write_trajectory_snapshot(self, fout=None):
        
        if fout is None:
            fout = sys.stdout

        print(f"# STMSNAP {self.STMStep}", file=fout)
        self.print_path_summary_data(fout)
        print("",file=fout) # necessary for gnuplot

# -------------------------------------------------------------------------

    def init_stm_statistics(self, args=None):
        """
        Initialise STM optimisation statistics and moving-average buffers.
        """

        self.CurrentPathLength = 0.0
        self.OldCPathLength = 0.0

        for bead in self.beads:
            bead.OPos[:] = bead.Pos[:]

        self.PLenChange = 0.0

        self.MaxBeadMove = 0.0
        self.MaxBeadMoveID = 0
        self.AveBeadMove = 0.0

        self.MaxpMFSize = 0.0
        self.MaxpMFSizeID = 0
        self.AvepMFSize = 0.0

        self.MAPLenChange = 0.0
        self.MAMaxBeadMove = 0.0
        self.MAAveBeadMove = 0.0
        self.MAMaxpMFSize = 0.0
        self.MAAvepMFSize = 0.0

        self.TermCrit = 0

        # Moving-average buffer length.
        self.MABufLength = args.mabuflen

        self.MABufPLenChange    = np.zeros(self.MABufLength, dtype=float)
        self.MABufMaxBeadMove   = np.zeros(self.MABufLength, dtype=float)
        self.MABufAveBeadMove   = np.zeros(self.MABufLength, dtype=float)
        self.MABufMaxpMFSize    = np.zeros(self.MABufLength, dtype=float)
        self.MABufAvepMFSize    = np.zeros(self.MABufLength, dtype=float)

        # Final convergence thresholds.
        self.FinalPLenChange    = args.final_plenchange
        self.FinalMaxBeadMove   = args.final_maxbeadmove
        self.FinalAveBeadMove   = args.final_avebeadmove
        self.FinalMaxpMFSize    = args.final_maxpmfsize
        self.FinalAvepMFSize    = args.final_avepmfsize

# -------------------------------------------------------------------------

    def print_stm_header_f(self, fout=None):
        """
        Print STM optimisation table header.
        """

        import sys

        if fout is None:
            fout = sys.stdout

        print("#", file=fout)
        print("# NOTE: All values are in scaled units.", file=fout)
        print("#", file=fout)
        print(
            "#     |                                         Current Values                                          |    |                             Running Averages                             |",
            file=fout,
        )
        print(
            "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|",
            file=fout,
        )
        print(
            "# Step|  Path length  Length change | Max bead move BID Ave bead move | Max pMF size  BID  Ave pMF size |Term|Length change | Max bead move Ave bead move | Max pMF size   Ave pMF size |",
            file=fout,
        )
        print(
            "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|",
            file=fout,
        )
        print(
            "#    1|             2              3|             4   5              6|             7   8              9|  10|            11|            12             13|            14             15|",
            file=fout,
        )
        print(
            "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|",
            file=fout,
        )

# -------------------------------------------------------------------------

    def calculate_stm_step_stat(self):
        """
        Calculate STM optimisation statistics.
        """

        self.PLenChange = self.CurrentPathLength - self.OldCPathLength

        self.MaxBeadMove = 0.0
        self.MaxBeadMoveID = 0
        self.AveBeadMove = 0.0

        self.MaxpMFSize = 0.0
        self.MaxpMFSizeID = 0
        self.AvepMFSize = 0.0

        bn = 0
        b = 0
        for bead in self.beads:
            b += 1
            if (self.freeterminals == False) and ( (b == 1) or (b == self.nbeads) ):
                continue
            
            # print(bead.Pos,bead.OPos)
            bmov = float(np.linalg.norm(bead.Pos - bead.OPos))
            mfsize = float(np.linalg.norm(bead.pGrad))

            self.AveBeadMove += bmov
            if bmov > self.MaxBeadMove:
                self.MaxBeadMove = bmov
                self.MaxBeadMoveID = b + 1

            self.AvepMFSize += mfsize
            if mfsize > self.MaxpMFSize:
                self.MaxpMFSize = mfsize
                self.MaxpMFSizeID = b + 1

            bn += 1

        if bn > 0:
            self.AveBeadMove /= float(bn)
            self.AvepMFSize /= float(bn)

        # Shift moving-average buffers to the left.
        self.MABufPLenChange[:-1] = self.MABufPLenChange[1:]
        self.MABufMaxBeadMove[:-1] = self.MABufMaxBeadMove[1:]
        self.MABufAveBeadMove[:-1] = self.MABufAveBeadMove[1:]
        self.MABufMaxpMFSize[:-1] = self.MABufMaxpMFSize[1:]
        self.MABufAvepMFSize[:-1] = self.MABufAvepMFSize[1:]

        # Add new values.
        self.MABufPLenChange[-1] = self.PLenChange
        self.MABufMaxBeadMove[-1] = self.MaxBeadMove
        self.MABufAveBeadMove[-1] = self.AveBeadMove
        self.MABufMaxpMFSize[-1] = self.MaxpMFSize
        self.MABufAvepMFSize[-1] = self.AvepMFSize

        # Reset moving averages.
        self.MAPLenChange = 0.0
        self.MAMaxBeadMove = 0.0
        self.MAAveBeadMove = 0.0
        self.MAMaxpMFSize = 0.0
        self.MAAvepMFSize = 0.0
        self.TermCrit = 0

        if self.STMStep < self.MABufLength:
            return

        self.MAPLenChange = float(np.mean(self.MABufPLenChange))
        self.MAMaxBeadMove = float(np.mean(self.MABufMaxBeadMove))
        self.MAAveBeadMove = float(np.mean(self.MABufAveBeadMove))
        self.MAMaxpMFSize = float(np.mean(self.MABufMaxpMFSize))
        self.MAAvepMFSize = float(np.mean(self.MABufAvepMFSize))

        # Determine number of fulfilled termination criteria.
        if self.MAPLenChange < self.FinalPLenChange:
            self.TermCrit += 1

        if self.MAMaxBeadMove < self.FinalMaxBeadMove:
            self.TermCrit += 1

        if self.MAAveBeadMove < self.FinalAveBeadMove:
            self.TermCrit += 1

        if self.MAMaxpMFSize < self.FinalMaxpMFSize:
            self.TermCrit += 1

        if self.MAAvepMFSize < self.FinalAvepMFSize:
            self.TermCrit += 1

# -------------------------------------------------------------------------

    def print_stm_step_info_f(self, fout=None):
        """
        Print one STM optimisation-statistics line.
        """

        import sys

        if fout is None:
            fout = sys.stdout

        print(
            f"{self.STMStep:>6d} "
            f"{self.CurrentPathLength:14.7e} "
            f"{self.PLenChange:14.7e} "
            f"{self.MaxBeadMove:14.7e} "
            f"{self.MaxBeadMoveID:3d} "
            f"{self.AveBeadMove:14.7e} "
            f"{self.MaxpMFSize:14.7e} "
            f"{self.MaxpMFSizeID:3d} "
            f"{self.AvepMFSize:14.7e} "
            f" {self.TermCrit:1d}/5 "
            f"{self.MAPLenChange:14.7e} "
            f"{self.MAMaxBeadMove:14.7e} "
            f"{self.MAAveBeadMove:14.7e} "
            f"{self.MAMaxpMFSize:14.7e} "
            f"{self.MAAvepMFSize:14.7e} ",
            file=fout,
        )

# =============================================================================
# CLI
# =============================================================================

def parse_args():

    def parse_csv_floats(text: str) -> np.ndarray:
        values = [float(v) for v in text.replace(";", ",").split(",") if v.strip()]
        if not values:
            raise argparse.ArgumentTypeError("at least one value is required")
        return np.asarray(values, dtype=float)

    parser = argparse.ArgumentParser(
            description="String-method path optimisation on a 2D FES/PES represented by Gaussian RBFs."
        )

    # -------------------------------------------------------------------------
    # CV1
    # -------------------------------------------------------------------------

    cv1group = parser.add_argument_group("The first collective variable (CV1) specification")

    cv1group.add_argument(
        "--cv1label", type=str, default=r"cv1",
        help="Label of the first collective variable."
    )

    cv1group.add_argument(
        "--cv1name", type=str, default=r"cv1",
        help="Name of the first collective variable."
    )

    cv1group.add_argument(
        "--cv1type", type=str, default=r"DIS",
        help="Type of the first collective variable."
    )

    cv1group.add_argument(
        "--cv1min", type=float, required=True,
        help="Minimum value of the first collective variable."
    )

    cv1group.add_argument(
        "--cv1max", type=float, required=True,
        help="Maximum value of the first collective variable."
    )

    cv1group.add_argument(
        "--cv1maxmove", type=float,
        help="Maximum move allowed for the first collective variable during the path optimization."
    )

    cv1group.add_argument(
        "--cv1nbins", type=int, required=True,
        help="Number of grid bins/points for the first collective variable."
    )

    # -------------------------------------------------------------------------
    # CV2
    # -------------------------------------------------------------------------

    cv2group = parser.add_argument_group("The second collective variable (CV2) specification")

    cv2group.add_argument(
        "--cv2label", type=str, default=r"cv2",
        help="Label of the second collective variable."
    )

    cv2group.add_argument(
        "--cv2name", type=str, default=r"cv2",
        help="Name of the second collective variable."
    )

    cv2group.add_argument(
        "--cv2type", type=str, default=r"DIS",
        help="Type of the second collective variable."
    )

    cv2group.add_argument(
        "--cv2min", type=float, required=True,
        help="Minimum value of the second collective variable."
    )

    cv2group.add_argument(
        "--cv2max", type=float, required=True,
        help="Maximum value of the second collective variable."
    )

    cv1group.add_argument(
        "--cv2maxmove", type=float,
        help="Maximum move allowed for the second collective variable during the path optimization."
    )

    cv2group.add_argument(
        "--cv2nbins", type=int, required=True,
        help="Number of grid bins/points for the second collective variable."
    )

    # -------------------------------------------------------------------------
    # Energy label
    # -------------------------------------------------------------------------

    enegroup = parser.add_argument_group("The energy axis specification")

    enegroup.add_argument(
        "--enen", type=str, default=r"${\Delta}G [kcal/mol]$",
        help="Energy label."
    )

    enegroup.add_argument(
        "--zmax", type=float, default=18.0,
        help="Maximum energy value considered."
    )

    enegroup.add_argument(
        "--contour_spacing", type=float, default=1.0,
        help="Contour spacing."
    )

    # -------------------------------------------------------------------------
    # RBF interpolation
    # -------------------------------------------------------------------------

    rbfgroup = parser.add_argument_group("The RBF (Radial Basis Function) interpolation specification")

    rbfgroup.add_argument(
        "--cv1nrbfs", type=int, default=20,
        help="Number of RBFs for the first collective variable."
    )

    rbfgroup.add_argument(
        "--cv2nrbfs", type=int, default=20,
        help="Number of RBFs for the second collective variable."
    )

    rbfgroup.add_argument(
        "--sx", type=float, default=1.5,
        help="Width factor for CV1 in the RBF static width mode."
    )

    rbfgroup.add_argument(
        "--sy", type=float, default=1.5,
        help="Width factor for CV2 in the RBF static width mode."
    )

    rbfgroup.add_argument(
        "--rcond", type=float, default=1.0e-9,
        help="SVD cutoff for RBF fitting.")

    # -------------------------------------------------------------------------
    # Files
    # -------------------------------------------------------------------------

    filegroup = parser.add_argument_group("The input/output files specification")

    filegroup.add_argument(
        "--input", required=True,
        help="Input FES/PES file. First three columns are CV1, CV2, energy."
    )

    filegroup.add_argument(
        "--output", default="_stm.path", 
        help="Path in the PMFLib format, printed at the beggining and end of STM."
    )
    
    filegroup.add_argument(
        "--summary", default="_stm.results", 
        help="Path summary, printed at the beggining and end of STM."
    )

    filegroup.add_argument(
        "--trajectory", default="_stm.traj", 
        help="Path summary trajectory."
    )

    filegroup.add_argument(
        "--optlog", default="_stm.log", 
        help="STM optimization log file."
    )
 
    # -------------------------------------------------------------------------
    # Path
    # -------------------------------------------------------------------------
    
    pathgroup = parser.add_argument_group("Path specification")

    pathgroup.add_argument("--initial-cv1", type=parse_csv_floats, required=True,
        help="Comma-separated initial CV1 path points."
    )

    pathgroup.add_argument("--initial-cv2", type=parse_csv_floats, required=True,
        help="Comma-separated initial CV2 path points."
    )

    pathgroup.add_argument("--pathname", type=str, default="p1",
        help="Name of the path."
    )


    pathgroup.add_argument("--nbeads", type=int, default=51,
        help="Number of string beads."
    )

    pathgroup.add_argument("--freeterminals", action="store_true", default=False,
        help="Optimise endpoints by steepest descent."
    )

    pathgroup.add_argument("--segdisc", type=int, default=10, 
        help="Path segment discrimination."
    )

    pathgroup.add_argument("--cvspline", type=int, default=1,
        help="Type of CV spline: 0 - interpolating cubic spline, 1 - smoothing cubic spline."
    )

    pathgroup.add_argument("--spline-lambda", type=float, default=0.999,
        help="Lambda for the internal smoothing cubic spline; 1.0 gives interpolation."
    )
    
    pathgroup.add_argument("--spline-sigma", type=float, default=0.01,
        help="Default sigma assigned to all knots of the internal smoothing cubic spline."
    )

    # -------------------------------------------------------------------------
    # STM Setup
    # -------------------------------------------------------------------------

    stmgroup = parser.add_argument_group("Path specification")

    stmgroup.add_argument("--nstepmax", type=int, default=200, 
        help="Maximum optimisation steps."
    )
    
    stmgroup.add_argument("--smoothingfac", type=float, default=0.0,
        help="Path smoothing factor."
    )
    
    stmgroup.add_argument("--smoothinterval", type=int, default=0,
        help="How often to smooth the path."
    )
    
    stmgroup.add_argument("--reparaminterval", type=int, default=1,
        help="How often to reparametrize the path."
    )

    # -------------------------------------------------------------------------
    # STM Setup
    # -------------------------------------------------------------------------

    adagroup = parser.add_argument_group("ADABelif specification")

    adagroup.add_argument("--stepsize", type=float, default=0.003,
        help="Optimisation time step."
    )

    adagroup.add_argument("--beta1", type=float, default=0.7,
        help="ADABelif beta1."
    )

    adagroup.add_argument("--beta2", type=float, default=0.99,
        help="ADABelif beta2."
    )

    adagroup.add_argument("--mingnormesp", type=float, default=1e-7,
        help="ADABelif epsilon."
    )

    # -------------------------------------------------------------------------
    # Termination
    # -------------------------------------------------------------------------

    termgroup = parser.add_argument_group("Termination criteria for the STM path optimization")

    termgroup.add_argument("--mabuflen", type=int, default=3,
        help="Moving-average buffer length."
    )

    termgroup.add_argument("--final-plenchange", type=float, default=0.001,
        help="Final threshold for moving-average path-length change."
    )
    
    termgroup.add_argument("--final-maxbeadmove", type=float, default=0.005,
        help="Final threshold for moving-average maximum bead movement."
    )

    termgroup.add_argument("--final-avebeadmove", type=float, default=.005,
        help="Final threshold for moving-average average bead movement."
    )
    
    termgroup.add_argument("--final-maxpmfsize", type=float, default=2.00,
        help="Final threshold for moving-average maximum projected mean-force size."
    )

    termgroup.add_argument("--final-avepmfsize", type=float, default=0.80,
        help="Final threshold for moving-average average projected mean-force size."
    )

    # -------------------------------------------------------------------------
    # Plots
    # -------------------------------------------------------------------------

    plotgroup = parser.add_argument_group("The graphical plot specification")

    plotgroup.add_argument("--plot", action="store_true", default=False,
        help="Write PNG plots."
        )
    
    plotgroup.add_argument("--plot-prefix", type=str, default="_stm", 
        help="Prefix for output PNG plots."
        )
    
    plotgroup.add_argument("--show", action="store_true", default=False,
        help="Show plots interactively after saving."
    )

    plotgroup.add_argument(
        "--dpi", type=int, default=300,
        help="Resolution for plot figures."
    )

    return parser.parse_args()

# ------------------------------------------------------------------------------

def main() -> None:

    print("#")
    print("# ==============================================================================")
    print("#              *** Simplified String Method on 2D Energy Surface ***            ")
    print("# ==============================================================================")
    print("#         The stm-path-2D-surface utility is the part of PMFLib toolkit.        ")
    print("#")
    print("#==============================================================================#")
    print("# PMFLib - Potential of Mean Force Toolkit                                     #")
    print("# -----------------------------------------------------------------------------#")
    print("# Authors: (c) 2026 Petr Kulhanek (NCBR)                                       #")
    print("#                                                                              #")
    print("# NCBR:    National Centre for Biomolecular Research, Masaryk University, CZ   #")
    print("#                                                                              #")
    print("# PMFLib is licensed under Lesser GPL v2.1 and above.                          #")
    print("#==============================================================================#")

    args = parse_args()

    options = vars(args)
    print(f"")
    print("# All arguments ...")
    print(options)

    # do all STM stuff :-)
    stmpath = STMPath(args)
    stmpath.stm_optimize(args)
    print("")
    
# -------------------------------------------------------------------------

if __name__ == "__main__":
    main()

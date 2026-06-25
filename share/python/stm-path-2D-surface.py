#!/usr/bin/env python3
"""
Simplified String Method on 2D Energy Surface

This is a Python rewrite of stm_path01.m and C++ implementation of STM in PMFLib.

The FES is represented by the same basic Gaussian RBF idea used in analyse-2D-surface.py: coordinates are scaled to
[0, 1], a tensor-product Gaussian RBF basis is fitted by SVD, and energies plus analytical gradients are evaluated from the fitted surface.

In this implementation, the metric tensor (MT) is considered as an unit matrix.
"""
# ==============================================================================

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

# ==============================================================================
# Cubic splines for path parametrisation
# ==============================================================================

class CVInterpolatingCubicSpline:
    """
    Natural interpolating cubic spline.

    The spline is represented on each interval [x_i, x_{i+1}] as

        S_i(t) = sd_i + sc_i*t + sb_i*t^2 + sa_i*t^3

    where

        t = alpha - x_i

    Public interface:
        update_points(alphas, values)
        spline(alpha)       -> value
        spline(alpha, 1)    -> first derivative
    """

    def __init__(self):
        self.clear()

# ------------------------------------------------------------------------------

    def clear(self):
        self.n = 0              # number of spline intervals
        self.x = np.zeros(0, dtype=float)
        self.y = np.zeros(0, dtype=float)

        # Polynomial coefficients per interval.
        self.sa = np.zeros(0, dtype=float)   # cubic coefficient
        self.sb = np.zeros(0, dtype=float)   # quadratic coefficient
        self.sc = np.zeros(0, dtype=float)   # linear coefficient
        self.sd = np.zeros(0, dtype=float)   # constant coefficient

# ------------------------------------------------------------------------------

    def allocate(self, numofknots: int):
        self.clear()

        numofknots = int(numofknots)
        if numofknots < 2:
            raise RuntimeError("at least two knots are required")

        self.n = numofknots - 1

        self.x = np.zeros(numofknots, dtype=float)
        self.y = np.zeros(numofknots, dtype=float)

        self.sa = np.zeros(self.n, dtype=float)
        self.sb = np.zeros(self.n, dtype=float)
        self.sc = np.zeros(self.n, dtype=float)
        self.sd = np.zeros(self.n, dtype=float)

# ------------------------------------------------------------------------------

    def set_point(self, knotid: int, alpha: float, cv: float):
        if self.x.size == 0:
            raise RuntimeError("spline storage is not allocated")

        if knotid < 0 or knotid >= self.x.size:
            raise RuntimeError("knotid is out-of-range")

        self.x[knotid] = float(alpha)
        self.y[knotid] = float(cv)

# ------------------------------------------------------------------------------

    def update_points(self, alphas, values):
        alphas = np.asarray(alphas, dtype=float)
        values = np.asarray(values, dtype=float)

        if alphas.ndim != 1 or values.ndim != 1:
            raise ValueError("alphas and values must be one-dimensional arrays")

        if alphas.size != values.size:
            raise ValueError("alphas and values must have the same length")

        if alphas.size < 2:
            raise RuntimeError("at least two knots are required")

        if np.any(np.diff(alphas) <= 0.0):
            raise RuntimeError("spline knots must have strictly increasing alpha values")

        if self.x.size != alphas.size:
            self.allocate(alphas.size)

        self.x[:] = alphas
        self.y[:] = values

        self.build_spline()

# ------------------------------------------------------------------------------

    def build_spline(self):
        """
        Build a natural cubic spline.

        Natural boundary conditions are used:

            S''(x_0) = 0
            S''(x_n) = 0
        """

        if self.n <= 0:
            raise RuntimeError("not enough knots")

        self.sa[:] = 0.0
        self.sb[:] = 0.0
        self.sc[:] = 0.0
        self.sd[:] = 0.0

        h = np.diff(self.x)

        if np.any(h <= 0.0):
            raise RuntimeError("spline knots must have strictly increasing alpha values")

        # Special case: only one interval -> straight line.
        if self.n == 1:
            self.sd[0] = self.y[0]
            self.sc[0] = (self.y[1] - self.y[0]) / h[0]
            self.sb[0] = 0.0
            self.sa[0] = 0.0
            return

        # Solve for second-derivative-related coefficients.
        #
        # We solve the tridiagonal system for c_i where
        #
        #     S_i(t) = a_i + b_i t + c_i t^2 + d_i t^3
        #
        # with natural boundary conditions c_0 = c_n = 0.
        #
        # The internal equations are:
        #
        #     h_{i-1} c_{i-1}
        #   + 2 (h_{i-1} + h_i) c_i
        #   + h_i c_{i+1}
        #   = 3 [ (y_{i+1}-y_i)/h_i - (y_i-y_{i-1})/h_{i-1} ]
        #
        # for i = 1, ..., n-1.

        lower = np.zeros(self.n - 1, dtype=float)
        diag = np.zeros(self.n - 1, dtype=float)
        upper = np.zeros(self.n - 1, dtype=float)
        rhs = np.zeros(self.n - 1, dtype=float)

        for i in range(1, self.n):
            k = i - 1

            lower[k] = h[i - 1] if k > 0 else 0.0
            diag[k] = 2.0 * (h[i - 1] + h[i])
            upper[k] = h[i] if k < self.n - 2 else 0.0

            rhs[k] = 3.0 * (
                (self.y[i + 1] - self.y[i]) / h[i]
                - (self.y[i] - self.y[i - 1]) / h[i - 1]
            )

        c_internal = self._solve_tridiagonal(lower, diag, upper, rhs)

        c = np.zeros(self.n + 1, dtype=float)
        c[1:self.n] = c_internal

        # Convert to the coefficient convention used by this class:
        #
        #     S_i(t) = sd_i + sc_i*t + sb_i*t^2 + sa_i*t^3

        for i in range(self.n):
            self.sd[i] = self.y[i]
            self.sc[i] = (
                (self.y[i + 1] - self.y[i]) / h[i]
                - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0
            )
            self.sb[i] = c[i]
            self.sa[i] = (c[i + 1] - c[i]) / (3.0 * h[i])

# ------------------------------------------------------------------------------

    @staticmethod
    def _solve_tridiagonal(lower, diag, upper, rhs):
        """
        Solve a tridiagonal linear system using the Thomas algorithm.

        lower[i] is the subdiagonal element in row i.
        diag[i]  is the diagonal element in row i.
        upper[i] is the superdiagonal element in row i.
        """

        n = rhs.size

        if n == 0:
            return np.zeros(0, dtype=float)

        a = lower.copy()
        b = diag.copy()
        c = upper.copy()
        d = rhs.copy()

        for i in range(1, n):
            if b[i - 1] == 0.0:
                raise RuntimeError("singular tridiagonal system")

            w = a[i] / b[i - 1]
            b[i] -= w * c[i - 1]
            d[i] -= w * d[i - 1]

        if b[-1] == 0.0:
            raise RuntimeError("singular tridiagonal system")

        x = np.zeros(n, dtype=float)
        x[-1] = d[-1] / b[-1]

        for i in range(n - 2, -1, -1):
            if b[i] == 0.0:
                raise RuntimeError("singular tridiagonal system")

            x[i] = (d[i] - c[i] * x[i + 1]) / b[i]

        return x

# ------------------------------------------------------------------------------

    def _interval_index(self, alpha: float) -> tuple[int, float]:
        if self.n <= 0:
            raise RuntimeError("not enough knots")

        a = float(np.clip(alpha, self.x[0], self.x[-1]))

        i = int(np.searchsorted(self.x, a, side="right") - 1)
        i = min(max(i, 0), self.n - 1)

        dx = a - self.x[i]

        return i, dx

# ------------------------------------------------------------------------------

    def get_cv(self, alpha: float) -> float:
        i, dx = self._interval_index(alpha)

        return float(
            self.sd[i]
            + self.sc[i] * dx
            + self.sb[i] * dx * dx
            + self.sa[i] * dx * dx * dx
        )

# ------------------------------------------------------------------------------

    def get_cv_first_der(self, alpha: float) -> float:
        i, dx = self._interval_index(alpha)

        return float(
            self.sc[i]
            + 2.0 * self.sb[i] * dx
            + 3.0 * self.sa[i] * dx * dx
        )

# ------------------------------------------------------------------------------

    def __call__(self, alpha, der: int = 0):
        if der not in (0, 1):
            raise ValueError("only value (der=0) and first derivative (der=1) are supported")

        arr = np.asarray(alpha, dtype=float)

        if arr.ndim == 0:
            if der == 0:
                return self.get_cv(float(arr))
            else:
                return self.get_cv_first_der(float(arr))

        if der == 0:
            return np.array([self.get_cv(a) for a in arr], dtype=float)
        else:
            return np.array([self.get_cv_first_der(a) for a in arr], dtype=float)

# ==============================================================================
# CVSmoothingCubicSpline
# ==============================================================================

class CVSmoothingCubicSpline:
    """
    Natural smoothing cubic spline without SVD.

    The spline minimizes

        E(z) = sum_i w_i * (y_i - z_i)^2
             + lam * integral (S''(x))^2 dx

    where

        w_i = 1 / N

    Interface
    ---------
    spline = CVSmoothingCubicSpline(lam)

    spline.update_points(x, y)

    value = spline(x)
    der1  = spline(x, 1)
    der2  = spline(x, 2)
    der3  = spline(x, 3)

    Notes
    -----
    lam = 0.0 gives a natural interpolating cubic spline.

    Larger lam gives stronger smoothing.

    This implementation does not use SVD. It assumes that the involved
    matrices are symmetric positive definite and solves them by Cholesky
    factorization when possible.
    """

    def __init__(self, lam: float):
        self.lam = float(lam)

        if self.lam < 0.0:
            raise ValueError("lam must be non-negative.")

        self.x = None
        self.y = None
        self.w = None

        self.npts = 0
        self.nseg = 0
        self.h = None

        # Smoothed knot values.
        self.z = None

        # Polynomial coefficients on interval i:
        #
        #   S_i(x) = a_i + b_i*t + c_i*t^2 + d_i*t^3
        #
        # where:
        #
        #   t = x - x_i
        #
        self.a = None
        self.b = None
        self.c = None
        self.d = None

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def update_points(self, x, y):
        """
        Set or update spline points.

        Parameters
        ----------
        x : array_like
            Knot positions. They must be strictly increasing after sorting.

        y : array_like
            Values at knot positions.

        Important
        ---------
        This method does not merge duplicate or nearly duplicate x values.
        Duplicate x values are treated as an error.
        """

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if x.ndim != 1:
            raise ValueError("x must be a one-dimensional array.")

        if y.ndim != 1:
            raise ValueError("y must be a one-dimensional array.")

        if x.size != y.size:
            raise ValueError("x and y must have the same length.")

        if x.size < 2:
            raise ValueError("At least two points are required.")

        order = np.argsort(x)
        x = x[order]
        y = y[order]

        h = np.diff(x)

        if np.any(h <= 0.0):
            raise ValueError("x values must be strictly increasing.")

        self.x = x
        self.y = y

        # Uniform weights:
        #
        #   w_i = 1 / N
        #
        self.w = np.full(x.size, 1.0 / float(x.size))

        self.npts = x.size
        self.nseg = x.size - 1
        self.h = h

        self._build_spline()

    def __call__(self, x_eval, der: int = 0):
        """
        Evaluate spline or derivative.

        Parameters
        ----------
        x_eval : float or array_like
            Evaluation point or points.

        der : int, default=0
            Derivative order.

            der = 0:
                spline value

            der = 1:
                first derivative

            der = 2:
                second derivative

            der = 3:
                third derivative

            der > 3:
                zero
        """

        if self.x is None:
            raise RuntimeError(
                "Spline is not initialized. Call update_points(x, y) first."
            )

        der = int(der)

        if der < 0:
            raise ValueError("Derivative order must be non-negative.")

        scalar_input = np.isscalar(x_eval)
        x_eval = np.asarray(x_eval, dtype=float)

        idx = np.searchsorted(self.x, x_eval, side="right") - 1
        idx = np.clip(idx, 0, self.nseg - 1)

        t = x_eval - self.x[idx]

        a = self.a[idx]
        b = self.b[idx]
        c = self.c[idx]
        d = self.d[idx]

        if der == 0:
            out = a + b*t + c*t**2 + d*t**3
        elif der == 1:
            out = b + 2.0*c*t + 3.0*d*t**2
        elif der == 2:
            out = 2.0*c + 6.0*d*t
        elif der == 3:
            out = 6.0*d
        else:
            out = np.zeros_like(x_eval, dtype=float)

        if scalar_input:
            return float(out)

        return out

    # -------------------------------------------------------------------------
    # Spline construction
    # -------------------------------------------------------------------------

    def _build_spline(self):
        """
        Build the smoothing spline.
        """

        if self.npts == 2:
            self.z = self.y.copy()
            self._build_natural_cubic_from_values(self.z)
            return

        if self.lam == 0.0:
            self.z = self.y.copy()
        else:
            self.z = self._calculate_smoothed_values()

        self._build_natural_cubic_from_values(self.z)

    def _calculate_smoothed_values(self):
        """
        Calculate smoothed knot values z from

            (W + lam*K) z = W y

        without using SVD.
        """

        q, r = self._build_reinsch_matrices()

        # Compute
        #
        #   K = Q R^{-1} Q^T
        #
        # without explicitly inverting R.
        #
        # Solve
        #
        #   R X = Q^T
        #
        # and then use
        #
        #   K = Q X
        #
        xmat = self._solve_spd(r, q.T)
        k = q @ xmat

        lhs = np.diag(self.w) + self.lam * k
        rhs = self.w * self.y

        z = self._solve_spd(lhs, rhs)

        return z

    def _build_reinsch_matrices(self):
        """
        Build Q and R matrices for the natural cubic smoothing spline.

        The roughness penalty can be written as

            z^T K z

        with

            K = Q R^{-1} Q^T
        """

        n = self.npts
        h = self.h

        q = np.zeros((n, n - 2), dtype=float)

        for i in range(n - 2):
            q[i, i] = 1.0 / h[i]
            q[i + 1, i] = -1.0 / h[i] - 1.0 / h[i + 1]
            q[i + 2, i] = 1.0 / h[i + 1]

        r = np.zeros((n - 2, n - 2), dtype=float)

        for i in range(n - 2):
            r[i, i] = (h[i] + h[i + 1]) / 3.0

        for i in range(n - 3):
            r[i, i + 1] = h[i + 1] / 6.0
            r[i + 1, i] = h[i + 1] / 6.0

        return q, r

    def _build_natural_cubic_from_values(self, z):
        """
        Build natural cubic spline coefficients through values z.

        On interval [x_i, x_{i+1}]:

            S_i(x) = a_i + b_i*t + c_i*t^2 + d_i*t^3

        where:

            t = x - x_i
        """

        n = self.npts
        h = self.h

        if n == 2:
            self.a = np.array([z[0]], dtype=float)
            self.b = np.array([(z[1] - z[0]) / h[0]], dtype=float)
            self.c = np.array([0.0], dtype=float)
            self.d = np.array([0.0], dtype=float)
            return

        amat = np.zeros((n - 2, n - 2), dtype=float)
        rhs = np.zeros(n - 2, dtype=float)

        for i in range(1, n - 1):
            row = i - 1

            if i > 1:
                amat[row, row - 1] = h[i - 1]

            amat[row, row] = 2.0 * (h[i - 1] + h[i])

            if i < n - 2:
                amat[row, row + 1] = h[i]

            rhs[row] = 6.0 * (
                (z[i + 1] - z[i]) / h[i]
                - (z[i] - z[i - 1]) / h[i - 1]
            )

        m_inner = self._solve_spd(amat, rhs)

        m = np.zeros(n, dtype=float)
        m[1:-1] = m_inner

        self.a = z[:-1].copy()

        self.b = (
            (z[1:] - z[:-1]) / h
            - h * (2.0*m[:-1] + m[1:]) / 6.0
        )

        self.c = m[:-1] / 2.0

        self.d = (m[1:] - m[:-1]) / (6.0*h)

    # -------------------------------------------------------------------------
    # Linear algebra helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _solve_spd(a, b):
        """
        Solve

            a x = b

        for a symmetric positive-definite matrix a.

        The preferred path is Cholesky factorization. If Cholesky fails,
        np.linalg.solve is tried as a fallback. No SVD is used.
        """

        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)

        try:
            lmat = np.linalg.cholesky(a)

            # Solve
            #
            #   L q = b
            #
            q = np.linalg.solve(lmat, b)

            # Solve
            #
            #   L.T x = q
            #
            x = np.linalg.solve(lmat.T, q)

            return x

        except np.linalg.LinAlgError:
            # Fallback for cases where the matrix is numerically not recognized
            # as positive definite. This still does not use SVD.
            return np.linalg.solve(a, b)

# ==============================================================================
# RBF surface model, adapted from analyse-2D-surface.py
# ==============================================================================

# ==============================================================================
# Axis
# ==============================================================================

class Axis:
    """One collective variable axis, represented internally on the scaled interval [0, 1]."""

    def __init__(self, cvmin, cvmax, nrbfs, npts, label):
        
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
        self.label      = str(label)

        self.range = self.cvmax - self.cvmin
        self.width = 1.0 / self.nrbfs
        self.width_scale = 1.0

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

# ==============================================================================
# EnergySurface2D
# ==============================================================================

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

    def load(self, filename: str | Path, xcolumn=1, ycolumn=2, ecolumn=3) -> None:
        """Load text FES data. First three columns are CV1, CV2, energy."""

        x_values = []
        y_values = []
        e_values = []

        with open(filename, "r", encoding="utf-8") as fin:
            for line in fin:
                line = line.strip()

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                line = line.split("#", 1)[0].strip()

                if not line:
                    continue

                fields = line.split()

                if len(fields) < max(xcolumn,ycolumn,ecolumn):
                    continue

                x_values.append(float(fields[xcolumn-1]))
                y_values.append(float(fields[ycolumn-1]))
                e_values.append(float(fields[ecolumn-1]))

        if len(e_values) == 0:
            raise ValueError("No valid data points were loaded")

        self.x_data = np.asarray(x_values, dtype=float)
        self.y_data = np.asarray(y_values, dtype=float)
        self.e_data = np.asarray(e_values, dtype=float)

        zmin = np.nanmin(self.e_data)
        self.e_data = self.e_data - zmin

        if self.zmax == None:
            self.zmax = np.nanmax(self.e_data)

        # determine sampling spacing -------------

        # sort and remove exact duplicates
        x_data_unique = np.unique(self.x_data)

        # differences between consecutive unique sorted values
        dx = np.diff(x_data_unique)

        # minimal positive difference
        self.x_min_dx = np.min(dx[dx > 0]) if np.any(dx > 0) else None

        # sort and remove exact duplicates
        y_data_unique = np.unique(self.y_data)

        # differences between consecutive unique sorted values
        dy = np.diff(y_data_unique)

        # minimal positive difference
        self.y_min_dy = np.min(dy[dy > 0]) if np.any(dy > 0) else None  

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

        x_edges = np.linspace(
            self.x_axis.cvmin,
            self.x_axis.cvmax,
            self.x_axis.npts + 1,
        )
        y_edges = np.linspace(
            self.y_axis.cvmin,
            self.y_axis.cvmax,
            self.y_axis.npts + 1,
        )

        # centered at bins
        x_grid = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_grid = 0.5 * (y_edges[:-1] + y_edges[1:])

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

        # ------------------------------------------------------------
        # Default threshold:
        # a grid point is sampled if there is a data point roughly within
        # one grid-cell diagonal.
        # ------------------------------------------------------------

        if max_distance is None:
            max_distance = 1.5 * np.sqrt(self.x_min_dx**2 + self.y_min_dy**2)

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
                            path_types: np.ndarray = None,
                            filename: str = None, show: bool = None,
                            figsize = None, dpi: int = 300, legend: bool = False) -> None:

            cmap, norm = self.make_colormap()
            cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing
            levels = np.arange(0.0, cmax + self.contour_spacing, self.contour_spacing)

            bead_colors = {
                "normal":    "green",
                "free":      "orange",
                "permanent": "black",
            }

            fig, ax = plt.subplots(figsize=figsize)
            im = ax.pcolormesh(self.X, self.Y, self.Z, shading="auto", cmap=cmap, norm=norm)
            fig.colorbar(im, ax=ax, label=self.ene_label)
            ax.contour(self.X, self.Y, self.Z, levels=levels, colors="black", linewidths=0.5)

            if path_x is not None and path_y is not None:
                path_x = np.asarray(path_x, dtype=float)
                path_y = np.asarray(path_y, dtype=float)

                if path_x.shape != path_y.shape:
                    raise ValueError("path_x and path_y must have the same shape")

                # Draw the connecting string first.  Bead markers are drawn later so
                # that their colours are not hidden by the line.
                ax.plot(path_x, path_y, "-", color="white", linewidth=1.0, zorder=4)

                if path_types is None:
                    path_types = np.full(path_x.shape, "unknown", dtype=object)
                else:
                    path_types = np.asarray(path_types, dtype=object)
                    if path_types.shape != path_x.shape:
                        raise ValueError("path_types must have the same shape as path_x and path_y")

                # Plot one bead type at a time.  This gives a clean legend and keeps
                # the type-to-colour mapping local to the plotting routine.
                for bead_type in list(bead_colors.keys()):
                    mask = path_types == bead_type
                    if not np.any(mask):
                        continue
                    ax.scatter(
                        path_x[mask], path_y[mask],
                        s=10,
                        c=bead_colors[bead_type],
                        edgecolors="white",
                        linewidths=0.3,
                        label=bead_type,
                        zorder=5,
                    )

                if legend:
                    ax.legend(loc="best", frameon=True)

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

# ==============================================================================
# String-method utilities
# ==============================================================================

# ==============================================================================
# Bead
# ==============================================================================

class Bead:
    def __init__(self, stmpath: STMPath):

        self.ncvs       = stmpath.ncvs
        self.eval_uv    = stmpath.surface.eval_uv

        self.type       = None
        self.OPos       = np.zeros(self.ncvs)               # old position, scaled
        self.Pos        = np.zeros(self.ncvs)               # bead position, scaled

        self.Grad       = np.zeros(self.ncvs)             # ENE gradient
        self.pGrad      = np.zeros(self.ncvs)             # force acting perpendicularly to the path
        self.uGrad      = np.zeros(self.ncvs)             # gradient to move bead

        # per path segment
        self.dCVdAlpha  = np.zeros(self.ncvs)        
        self.Seglength  = 0.0

        self.Alpha      = None      # path position
        self.dAdAlpha   = None      # free energy derivative
        self.A          = None      # free energy, integrated
        self.Asurf      = None      # free energy from EnergySurface2D


# ------------------------------------------------------------------------------

    def UpdatePositionGradientDescent(self,stepsize,cvs):

        if self.type == "permanent":
            return
        
        for i in range(self.ncvs):

            dm = - stepsize * self.uGrad[i]

            if (cvs[i].smaxmov <= 0) or (math.fabs(dm) < cvs[i].smaxmov):
                self.Pos[i] = self.Pos[i] + dm
            else:
                self.Pos[i] = self.Pos[i] + cvs[i].smaxmov*math.copysign(1.0,dm)

# ==============================================================================
# STMPath
# ==============================================================================

class STMPath:
    def __init__(self, args):

        # path parameters
        self.ncvs           = 2
        # self.nbeads       -> load_path_file()
        # self.PathName     -> load_path_file()

        # axes/CVs
        self.cvs = []

        x_axis = Axis(args.cv1min,args.cv1max,args.cv1nrbfs,args.cv1nbins,args.cv1label)
        self.cvs.append(x_axis)

        y_axis = Axis(args.cv2min,args.cv2max,args.cv2nrbfs,args.cv2nbins,args.cv2label)
        self.cvs.append(y_axis)

        self.surface = EnergySurface2D(x_axis,y_axis,args.enelabel,args.zmax,args.contour_spacing)

        print("")
        print(f"# Load FES: {args.input_fes}")
        self.surface.load(args.input_fes,xcolumn=args.input_fes_x_column,ycolumn=args.input_fes_y_column,ecolumn=args.input_fes_e_column)

        print(f"")
        print(f"# Optimize RBF ...")
        print(f"  Sx:           {args.rbfsx:10.3f}")
        print(f"  Sy:           {args.rbfsy:10.3f}")

        self.surface.fit(sx=args.rbfsx, sy=args.rbfsy, rcond=args.rcond)
        print(f"  RBF fit RMSE: {self.surface.fit_rmse:10.3f}")

        print("")
        print("# Calculate Z surface ...")
        self.surface.calc_grid()

        print("")
        print("# CV splines ...")
        if args.cvspline == 1:
            print(f"  >>> Smoothing Cubic Spline")
            print(f"      Lambda: {args.spline_lambda:10.5e}")
            self.cv_splines =  [
                CVSmoothingCubicSpline(lam=args.spline_lambda)
                for _ in range(self.ncvs)
            ]
        else:
            print(f"  >>> Interpolating Cubic Spline")
            self.cv_splines =  [CVInterpolatingCubicSpline() for _ in range(self.ncvs)]

        print("")
        print("# Initial path ...")

        # read initial path
        input_beads = self.load_path_file(args.input_path)

        # plot the user initial path
        if args.plot :
            path_x, path_y, path_types = self.beads_to_xy_arrays(input_beads)
            self.surface.plot_fes_with_path(title="Initial Path - User Input",
                path_x=path_x, path_y=path_y, path_types=path_types,
                filename=f"{args.plot_prefix}_path_0000_a_initial-user.png", show=args.show, figsize=args.figsize, dpi=args.dpi)

        #  generate completed path
        self.beads = self.generate_beads_from_input_beads(input_beads)

        # plot the user completed path
        if args.plot :
            path_x, path_y, path_types = self.beads_to_xy_arrays(self.beads)
            self.surface.plot_fes_with_path(title="Initial Path - Full Path",
                path_x=path_x, path_y=path_y, path_types=path_types,
                filename=f"{args.plot_prefix}_path_0000_b_initial-full.png", show=args.show, figsize=args.figsize, dpi=args.dpi)

        # smoothing, reparameterization
        self.SmoothInterval         = args.smoothinterval
        self.SmoothingFac           = args.sfac
        self.ReparamInterval        = args.reparaminterval
        self.DetectMinimaAtStep     = args.detect_minima_at_step
        self.MinimaTreshold         = args.minima_threshold

        # STM
        self.STMStep            = 0
        self.StepSize           = args.stepsize
        self.nstepmax           = args.nstepmax


        # run initial statistics
        self.init_stm_statistics(args)
        
    # --------------------------------------------------------------------------

    def stm_optimize(self,args):

        print("")
        print("# STM path optimization ...")

        self.calc_beads()
        self.integrate_path()
        self.calculate_stm_step_stat()

        if args.output_path is not None:
            with open(args.output_path,"w") as fout:
                self.print_path(fout=fout,beads=self.beads)

        print("")
        print("# Initial path summary ...")
        self.print_path_summary_header()
        self.print_path_summary_data()

        if args.summary is not None:
            with open(args.summary,"w") as fsum:
                self.print_path_summary_header(fout=fsum)
                self.print_path_summary_data(fout=fsum)

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

            if self.STMStep == self.DetectMinimaAtStep:
                # mark minima as free beads
                self.detect_minima();
                self.StepSize *= args.scale_stepsize
                print(f"# >>>>> New Stepsize: {self.StepSize}")
                print(f"# >>>>> CV Splines ...")
                if args.cvspline == 1:
                    print(f"#       Smoothing Cubic Spline")
                    print(f"#       Lambda: {args.spline_lambda*args.scale_lambda:10.5e}")
                    self.cv_splines =  [
                        CVSmoothingCubicSpline(lam=args.spline_lambda*args.scale_lambda)
                        for _ in range(self.ncvs)
                    ]
                else:
                    print(f"  >>> Interpolating Cubic Spline")
                    self.cv_splines =  [CVInterpolatingCubicSpline() for _ in range(self.ncvs)]

            self.calculate_stm_step_stat()

            self.print_stm_step_info_f()
            if args.optlog is not None:
                self.print_stm_step_info_f(fout=flog)

            if args.trajectory is not None:
                self.write_trajectory_snapshot(ftraj)

            if args.plot :
                path_x, path_y, path_types = self.beads_to_xy_arrays(self.beads)
                self.surface.plot_fes_with_path(title=f"Intermediate Path #{self.STMStep:04d}",
                    path_x=path_x, path_y=path_y, path_types=path_types,
                    filename=f"{args.plot_prefix}_path_{self.STMStep:04d}_c.png", show=args.show, figsize=args.figsize, dpi=args.dpi)

            if self.TermCrit == 5:
                break

        if self.TermCrit == 5:
            print(">>> STM path optimization was successful.")
        else:
            print(">>> STM path optimization failed.")

        print("")
        print("# Final path summary ...")
        self.print_path_summary_header()
        self.print_path_summary_data()

        # now overwrite the summary with final data
        if args.summary is not None:
            with open(args.summary,"w") as fsum:
                self.print_path_summary_header(fout=fsum)
                self.print_path_summary_data(fout=fsum)

        if args.plot:
            path_x, path_y, path_types = self.beads_to_xy_arrays(self.beads)
            self.surface.plot_fes_with_path(title="Final Path",
                path_x=path_x, path_y=path_y, path_types=path_types,
                filename=f"{args.plot_prefix}_path_{self.STMStep:04d}_d_final.png", show=args.show, figsize=args.figsize, dpi=args.dpi)

        if args.output_path is not None:
            with open(args.output_path,"w") as fout:
                self.print_path(fout=fout,beads=self.beads)

        if args.optlog is not None:
            flog.close()

        if args.trajectory is not None:
            ftraj.close()

# ==============================================================================
# path manipulation
# ==============================================================================

    def parametrize_path(self, beads):
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

        total_length      = 0.0

        for b in range(1, len(beads)):
            diff = beads[b].Pos - beads[b - 1].Pos
            slen = np.linalg.norm(diff)

            total_length      += slen

        if total_length == 0.0:
            raise RuntimeError("path has zero length")

        # ---------------------------------------------------------------------
        # Initial alpha values from linear interpolation.
        # ---------------------------------------------------------------------

        beads[0].Alpha = 0.0

        path_length = 0.0

        for b in range(1, len(beads) - 1):
            diff = beads[b].Pos - beads[b - 1].Pos
            slen = np.linalg.norm(diff)

            if slen == 0.0:
                raise RuntimeError("path segment has zero length")
            
            path_length   += slen

            beads[b].Alpha = path_length / total_length

        beads[-1].Alpha = 1.0

        # ---------------------------------------------------------------------
        # Build the splines
        # ---------------------------------------------------------------------

        alphas = np.array([bead.Alpha for bead in beads], dtype=float)

        for i in range(self.ncvs):
            values = np.array([bead.Pos[i] for bead in beads], dtype=float)
            self.cv_splines[i].update_points(alphas, values)

        return total_length

    # --------------------------------------------------------------------------

    def update_all_positions(self, optimizer = 0):

        # backup old positions
        for bead in self.beads:
            bead.OPos[:]  = bead.Pos[:]

        self.OldCPathLength = self.CurrentPathLength

        self.STMStep += 1

        # update positions
        for bead in self.beads:
            bead.UpdatePositionGradientDescent(self.StepSize,self.cvs)

    # --------------------------------------------------------------------------

    def smooth_all_positions(self):

        if (self.SmoothInterval == 0) or (self.STMStep % self.SmoothInterval != 0 ):
            return;
    
        old_pos = np.array([bead.Pos.copy() for bead in self.beads])

        for i in range(1, self.nbeads - 1):
            if self.beads[i].type == "normal":
                self.beads[i].Pos[:] = (
                    (1.0 - self.SmoothingFac) * old_pos[i]
                    + 0.5 * self.SmoothingFac * (old_pos[i - 1] + old_pos[i + 1])
                )

    # --------------------------------------------------------------------------

    def iter_path_segments(self, beads):
        """
        Yield segments composed of:
            (left terminal bead) - optional, first path bead or free bead
            one or more normal / permanent beads
            (right terminal bead) - optional, last path bead or free bead
        """

        nbeads = len(beads)
        i = 0

        while i < nbeads:

            if i != 0 and beads[i].type != "free":
                raise RuntimeError(f"path segment must start with the free bead or the first bead of the path, bidx: {i}")

            seg_first = i

            i += 1

            # Find the end of the segment
            while i < nbeads and beads[i].type != "free":
                i += 1

            seg_last = i

            # make a list
            segment_beads = []
            for idx in range(seg_first,seg_last+1):
                if idx < 0 or idx >= self.nbeads:
                    continue
                segment_beads.append(self.beads[idx])

            yield segment_beads

    # --------------------------------------------------------------------------

    def reparametrize_all_positions(self):

        if (self.ReparamInterval == 0) or (self.STMStep % self.ReparamInterval != 0):
            return

        for segment_beads in self.iter_path_segments(self.beads):
            
            # At least one bead inside the segment
            if len(segment_beads) < 3:
                continue

            # Parametrize this local segment and build CV splines for it.
            # The first and last beads become alpha = 0.0 and alpha = 1.0.
            self.parametrize_path(segment_beads)

            # generate evenly distributed set of alphas
            alphas = []
            for idx, bead in enumerate(segment_beads):
                alpha = float(idx) / float(len(segment_beads) - 1)
                alphas.append(alpha)

            alphas[0] = 0.0
            alphas[-1] = 1.0

            # update positions along the path based on new alphas
            for cv in range(self.ncvs):
                # redistribute beads - exclude terminals
                for bidx, bead in enumerate(segment_beads[1:-1],start=1):
                    if bead.type == "normal":
                        bead.Pos[cv] = self.cv_splines[cv](alphas[bidx])
       
    # --------------------------------------------------------------------------

    def check_boundaries_of_beads(self, beads):
        """
        Correct bead.Pos positions after path generation.
        """
        for bead in beads:
            for cvidx, cv in enumerate(self.cvs):
                if bead.Pos[cvidx] < cv.scale(cv.pathmin):
                    bead.Pos[cvidx] = cv.scale(cv.pathmin)
                if bead.Pos[cvidx] > cv.scale(cv.pathmax):
                    bead.Pos[cvidx] = cv.scale(cv.pathmax)  

# ------------------------------------------------------------------------------

    def detect_minima(self):

        print(f"# >>>>> Detecting Minima ...")

        alphas = np.array([bead.Alpha for bead in self.beads], dtype=float)
        enes   = np.array([bead.Asurf for bead in self.beads], dtype=float)

        _, bidxs = self.find_minima_positions(alphas,enes,self.MinimaTreshold)

        for bidx in bidxs:
            self.beads[bidx].type = "free"

        print(f"#       Found: {len(bidxs)}")

# ------------------------------------------------------------------------------

    def find_minima_positions(self, x, y, ethr=0.0):
        """
        Detect minima of y = f(x), merge minima within the same basin,
        and return only the lowest minimum from each basin.

        Two neighbouring minima are considered separate only if the maximum
        between them is at least `ethr` higher than both minima. Otherwise,
        they belong to the same basin and only the lower one is retained.

        Parameters
        ----------
        x : array-like
            Positions.
        y : array-like
            Function values.
        ethr : float
            Minimum barrier height required to separate two basins.

        Returns
        -------
        minima_x : np.ndarray
            x positions of selected minima.
        minima_indices : list
            Indices or index ranges corresponding to selected minima.
        """

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("x and y must be one-dimensional arrays.")

        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")

        n = len(y)
        if n < 3:
            return np.array([]), []

        # ------------------------------------------------------------
        # Step 1: detect all local minima, including flat minima.
        # Each minimum is represented as (start, end, xmin, ymin).
        # For sharp minima, start == end.
        # ------------------------------------------------------------

        candidates = []

        i = 1
        while i < n - 1:

            # sharp minimum
            if y[i] < y[i - 1] and y[i] < y[i + 1]:
                candidates.append((i, i, x[i], y[i]))
                i += 1
                continue

            # flat minimum
            if y[i] < y[i - 1] and y[i] == y[i + 1]:
                start = i

                while i + 1 < n and y[i + 1] == y[start]:
                    i += 1

                end = i

                if end < n - 1 and y[end] < y[end + 1]:
                    xmin = 0.5 * (x[start] + x[end])
                    ymin = y[start]
                    candidates.append((start, end, xmin, ymin))

            i += 1

        if not candidates:
            return np.array([]), []

        # ------------------------------------------------------------
        # Step 2: merge minima that are not separated by ethr barrier.
        # ------------------------------------------------------------

        basins = []
        current_basin = [candidates[0]]

        for m1, m2 in zip(candidates[:-1], candidates[1:]):
            _, end1, _, ymin1 = m1
            start2, _, _, ymin2 = m2

            # maximum between the two minima
            ymax_between = np.max(y[end1:start2 + 1])

            barrier1 = ymax_between - ymin1
            barrier2 = ymax_between - ymin2

            separated = (barrier1 >= ethr) and (barrier2 >= ethr)

            if separated:
                basins.append(current_basin)
                current_basin = [m2]
            else:
                current_basin.append(m2)

        basins.append(current_basin)

        # ------------------------------------------------------------
        # Step 3: retain only the lowest minimum from each basin.
        # ------------------------------------------------------------

        selected = []

        for basin in basins:
            lowest = min(basin, key=lambda item: item[3])
            selected.append(lowest)

        minima_x = []
        minima_indices = []

        for start, end, xmin, ymin in selected:
            minima_x.append(xmin)

            if start == end:
                minima_indices.append(start)
            else:
                minima_indices.append((start, end))

        return np.asarray(minima_x), minima_indices
        
# ------------------------------------------------------------------------------

    def calc_beads(self):
        """
        For all beads:
            Calculate perpendicular projectors and projected mean forces
            Calculate kink angles
        
        Assumptions
        -----------
        - bead.Pos, bead.Grad, bead.pGrad are in scaled coordinates.
        - MTZ is the identity matrix and is omitted.
        - No additional CV range scaling is applied.
        - self.cv_splines were already built by parametrize_path()
        """

        # for all beads
        for bead in self.beads:            
            # Calculate MF and A
            bead.Asurf,grad,_ = self.surface.eval_uv(bead.Pos)
            bead.Grad[:] = grad[:]

        # for path segments
        for segment_beads in self.iter_path_segments(self.beads):
            
            # Parametrize this local segment and build CV splines for it.
            # The first and last beads become alpha = 0.0 and alpha = 1.0.
            seglen = self.parametrize_path(segment_beads)

            for bead in segment_beads:            

                bead.Seglength = seglen

                # Calculate tangent dCV/dalpha from the path splines.
                for i in range(self.ncvs):
                    bead.dCVdAlpha[i] = self.cv_splines[i](bead.Alpha, 1)

                slen2 = float(np.dot(bead.dCVdAlpha, bead.dCVdAlpha))
                if slen2 == 0.0:
                    raise RuntimeError("derivative segment has zero length")

                # get perpendicular gradient
                bead.pGrad[:] = bead.Grad[:] - bead.dCVdAlpha[:] * (np.dot(bead.dCVdAlpha, bead.Grad) / slen2)

                if bead.type == "free":
                    bead.uGrad[:] = bead.Grad[:]    # switch to gradient descent move
                elif bead.type == "normal":
                    bead.uGrad[:] = bead.pGrad[:]   # use perpendicular gradient
                else:
                    bead.uGrad[:] = 0.0

        # and now for the entire path
        self.CurrentPathLength = self.parametrize_path(self.beads)  

        # Calculate kink angles
        v1 = np.zeros(self.ncvs)
        v2 = np.zeros(self.ncvs)

        for bidx in range(self.nbeads):

            if bidx == 0  or bidx == self.nbeads - 1:
                self.beads[bidx].kangle = 180.0
                continue

            v1[:] = self.beads[bidx+1].Pos[:] - self.beads[bidx].Pos[:]
            v2[:] = self.beads[bidx-1].Pos[:] - self.beads[bidx].Pos[:]

            # 1. Compute dot product and magnitudes
            dot_product = np.dot(v1, v2)
            norm_v1 = np.linalg.norm(v1)
            norm_v2 = np.linalg.norm(v2)

            # 2. Get the cosine of the angle
            cos_theta = dot_product / (norm_v1 * norm_v2)

            # 3. Prevent floating-point errors from pushing cos_theta outside [-1, 1]
            cos_theta = np.clip(cos_theta, -1.0, 1.0)

            # 4. Calculate the angle in radians
            angle_radians = np.arccos(cos_theta)

            # 5. Convert to degrees (optional)
            self.beads[bidx].kangle = np.degrees(angle_radians)

# ------------------------------------------------------------------------------

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

            # bead.dCVdAlpha are per path segment, thus correction * self.CurrentPathLength / bead.Seglength
            bead.dAdAlpha = float(np.dot(bead.dCVdAlpha, bead.Grad)) * self.CurrentPathLength / bead.Seglength

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

        # shift global minimum to zero.
        amin = min(bead.A for bead in self.beads)
        emin = min(bead.Asurf for bead in self.beads)

        for bead in self.beads:
            bead.A -= amin
            bead.Asurf -= emin


   # --------------------------------------------------------------------------

    def generate_beads_from_input_beads(self, input_beads):
        """
        Generate the final path from user provided one.

        Parameters
        ----------
        input_beads : list[Bead]
            User-defined beads

        Returns
        -------
        beads : list[Bead]
            Final bead list.
        """

        # ---------------------------------------------------------------------
        # 1) Check input
        # ---------------------------------------------------------------------

        if len(input_beads) < 2:
            raise RuntimeError("At least two input beads (normal/free/permanent) must be specified in the input PATH file!")
        
        if self.nbeads < 2:
            raise RuntimeError("At least two beads (nbeads) must be requested in the input PATH file!")

        if len(input_beads) == self.nbeads:
            # check boundaries
            self.check_boundaries_of_beads(input_beads)

            # path is complete, keep it as it is
            return input_beads
        
        # path is incompleted, it will be rebuilded

        for bead in input_beads[1:-1]:
            if bead.type != "normal":
                raise RuntimeError(f"For the incomplete path, all inner beads must be 'normal' but {bead.type} was requested!")

        # ---------------------------------------------------------------------
        # 2) Parametrize user provided path
        # ---------------------------------------------------------------------

        self.parametrize_path(input_beads)

        # ---------------------------------------------------------------------
        # 3) Generate missing beads from the optimized input spline.
        # ---------------------------------------------------------------------

        beads = [Bead(self) for _ in range(self.nbeads)]

        for bidx, bead in enumerate(beads):
            alpha = float(bidx) / float(self.nbeads - 1)

            bead.Alpha = alpha
            bead.type  = "normal"

            # copy type of terminals from the input path
            if bidx == 0:
                bead.type = input_beads[0].type
            if bidx == self.nbeads - 1:
                bead.type = input_beads[-1].type

            for i in range(self.ncvs):
                bead.Pos[i] = self.cv_splines[i](alpha)

        # Force exact endpoint alphas.
        beads[0].Alpha = 0.0
        beads[-1].Alpha = 1.0

        # ---------------------------------------------------------------------
        # 4) Check boundaries.
        # ---------------------------------------------------------------------

        self.check_boundaries_of_beads(beads)

        return beads
    
    # --------------------------------------------------------------------------

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

        bead_types = np.array(
            ["unknown" if bead.type is None else str(bead.type) for bead in beads],
            dtype=object,
        )

        return xy[:, 0], xy[:, 1], bead_types

# ==============================================================================
# load/print path
# ==============================================================================

    def load_path_file(self,filename: str | Path):
        """
        Parse a PMFLib-like [PATH] file.

        Supported keywords
        ------------------
        [PATH]
        name
        nbeads
        ncvs
        names
        types
        min
        max
        maxmov
        normal      - moves only perpendicular to the path
        free        - moves by full gradient descent
        permanent   - does not move   
        """

    # --------------------------------------------------------------------------

        def strip_path_comment(line: str) -> str:
            """
            Remove PMFLib-style comments.

            Comment characters are '#', '!', and '*'.
            The first occurrence of any of them starts a comment.
            """

            cut = len(line)

            for c in ("#", "!", "*"):
                i = line.find(c)
                if i >= 0:
                    cut = min(cut, i)

            return line[:cut].strip()

    # --------------------------------------------------------------------------

        in_path_section = False

        name_loaded = False
        nbeads_loaded = False
        ncvs_loaded = False
        names_loaded = False
        types_loaded = False
        pathmin_loaded = False
        pathmax_loaded = False
        maxmov_loaded = False

        beads = []

        with open(filename, "r", encoding="utf-8") as fin:
            for lineno, raw_line in enumerate(fin, start=1):

                line = strip_path_comment(raw_line)

                if not line:
                    continue

                if line.upper() == "[PATH]":
                    in_path_section = True
                    continue

                if line.startswith("[") and line.endswith("]"):
                    in_path_section = False
                    continue

                if not in_path_section:
                    continue

                fields = line.split()
                if not fields:
                    continue

                key = fields[0].lower()
                values = fields[1:]

                try:
                    if key == "name":
                        if len(values) != 1:
                            raise ValueError("keyword 'name' expects one value")
                        self.PathName = values[0]
                        name_loaded = True

                    elif key == "nbeads":
                        if len(values) != 1:
                            raise ValueError("keyword 'nbeads' expects one value")
                        self.nbeads = int(values[0])
                        nbeads_loaded = True

                    elif key == "ncvs":
                        if len(values) != 1:
                            raise ValueError("keyword 'ncvs' expects one value")
                        if int(values[0]) != self.ncvs:
                            raise ValueError(f"exactly {self.ncvs} CVs must be specified in the path")
                        ncvs_loaded = True

                    elif key == "names":
                        if len(values) != self.ncvs:
                            raise ValueError(f"keyword 'names' expects {self.ncvs} values")
                        for idx, name in enumerate(values):
                            self.cvs[idx].name = name
                        names_loaded = True

                    elif key == "types":
                        if len(values) != self.ncvs:
                            raise ValueError(f"keyword 'types' expects {self.ncvs} values")
                        for idx, type in enumerate(values):
                            self.cvs[idx].type = type
                        types_loaded = True

                    elif key == "min":
                        if len(values) != self.ncvs:
                            raise ValueError(f"keyword 'min' expects {self.ncvs} values")
                        for idx, min in enumerate(values):
                            if float(min) < self.cvs[idx].cvmin:
                                raise ValueError(f"'min' value for CV {idx+1} must be within FES area: {float(min)} < {self.cvs[idx].cvmin}!")
                            if float(min) > self.cvs[idx].cvmax:
                                raise ValueError(f"p'min' value for CV {idx+1} must be within FES area: {float(min)} > {self.cvs[idx].cvmax}!")
                            self.cvs[idx].pathmin = float(min)
                        pathmin_loaded = True

                    elif key == "max":
                        if len(values) != self.ncvs:
                            raise ValueError(f"keyword 'max' expects {self.ncvs} values")
                        for idx, max in enumerate(values):
                            if float(min) < self.cvs[idx].cvmin:
                                raise ValueError(f"'max' value for CV {idx+1} must be within FES area: {float(max)} < {self.cvs[idx].cvmin}!")
                            if float(min) > self.cvs[idx].cvmax:
                                raise ValueError(f"p'max' value for CV {idx+1} must be within FES area: {float(max)} > {self.cvs[idx].cvmax}!")
                            self.cvs[idx].pathmax = float(max)
                        pathmax_loaded = True

                    elif key == "maxmov":
                        if len(values) != self.ncvs:
                            raise ValueError(f"keyword 'maxmov' expects {self.ncvs} values")
                        for idx, maxmov in enumerate(values):
                            self.cvs[idx].maxmov = float(maxmov)
                            self.cvs[idx].smaxmov = self.cvs[idx].maxmov / (self.cvs[idx].cvmax - self.cvs[idx].cvmin)
                        maxmov_loaded = True

                    elif key in ("permanent", "normal", "free"):
                        if  len(values) != self.ncvs:
                            raise ValueError(
                                f"keyword '{key}' expects {self.ncvs} CV values, "
                                f"got {len(values)}"
                            )

                        bead = Bead(self)
                        bead.type = key
                        for idx, pos in enumerate(values):
                            bead.Pos[idx] = self.cvs[idx].scale(float(pos))

                        beads.append(bead)
                    else:
                        raise ValueError(f"unknown keyword '{key}'")

                except ValueError as e:
                    raise ValueError(
                        f"{filename}:{lineno}: invalid [PATH] line: {raw_line.rstrip()}\n"
                        f"Reason: {e}"
                    ) from e

        if not (name_loaded and nbeads_loaded and ncvs_loaded and names_loaded and types_loaded and pathmin_loaded and pathmax_loaded and maxmov_loaded):
            raise ValueError(f"mandatory path item not loaded | name:{name_loaded}/nbeads:{nbeads_loaded}/ncvs:{ncvs_loaded}/names:{names_loaded}/types:{types_loaded}/min:{pathmin_loaded}/max:{pathmax_loaded}/maxmov:{maxmov_loaded}")
        
        if len(beads) < 2:
            raise ValueError(f"at least two beads must be provided | normal/free/permanent")
        
        return beads

    # --------------------------------------------------------------------------

    def print_path(self, fout=None, beads=None):
            """
            Print the STM path in the same format as CSTMPath::PrintPath().

            Notes
            -----
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
                print(f" {cv.pathmin:12.5e}", end="", file=fout)
            print(file=fout)

            print("max       ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.pathmax:12.5e}", end="", file=fout)
            print(file=fout)

            print("maxmov    ", end="", file=fout)
            for cv in self.cvs:
                print(f" {cv.maxmov:12.5e}", end="", file=fout)
            print(file=fout)

            for index, bead in enumerate(beads):
                print(f"{bead.type:<9} ", end="", file=fout)

                for i, cv in enumerate(self.cvs):
                    scaled = bead.Pos[i]
                    unscaled = cv.unscale(scaled)
                    print(f" {float(unscaled):12.5e}", end="", file=fout)

                print(file=fout)

# ==============================================================================
# path summary and summary trajectory
# ==============================================================================

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
        print("#  ID   Type  MO ST KinkA  alpha    dA/dalpha        Asurf     CID Updates", end="", file=fout)
        for i in range(self.ncvs):
            print(f"          CV{i + 1:<1d}", end="", file=fout)
        for i in range(self.ncvs):
            print(f"         sCV{i + 1:<1d}", end="", file=fout)
        for i in range(self.ncvs):
            print(f"     dA/dsCV{i + 1:<1d}", end="", file=fout)
        for i in range(self.ncvs):
            print(f" dsCV{i + 1:<1d}/dalpha", end="", file=fout)
        for i in range(self.ncvs):
            print(f"  -|F{i + 1:<1d}/dalpha", end="", file=fout)

        print(f"         Aint", end="", file=fout)
        print(file=fout)

        # Delimiters.
        delimiter = "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------"
        delimiter += " ------------" * (5 * self.ncvs)
        delimiter += " ------------"
        print(delimiter, file=fout)

        # CV metadata rows.
        print(f"{'#      names':<74}", end="", file=fout)
        for _ in range(5):
            for cv in self.cvs:
                print(f" {cv.name:>12}", end="", file=fout)
        print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      types':<74}", end="", file=fout)
        for _ in range(2):
            for cv in self.cvs:
                print(f" {cv.type:>12}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      min':<74}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.pathmin:12.5e}", end="", file=fout)
        for cv in self.cvs:
            print(f" {0.0:12.5e}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      max':<74}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.pathmax:12.5e}", end="", file=fout)
        for cv in self.cvs:
            print(f" {1.0:12.5e}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      maxmov':<74}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.maxmov:12.5e}", end="", file=fout)
        for _ in range(4):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(f"             ", end="", file=fout)
        print(file=fout)

        print(delimiter, file=fout)

        delimiter2 = "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------"
        delimiter2 += " uuuuuuuuuuuu" * (1 * self.ncvs)
        delimiter2 += " ssssssssssss" * (4 * self.ncvs)
        delimiter2 += " ------------"

        print(delimiter2, file=fout)
        print(delimiter, file=fout)

        print("#    1      2  3  4     5      6            7            8       9      10", end="", file=fout)
        column_id = 11
        for _ in range(5):
            for _ in range(self.ncvs):
                print(f"{column_id:13d}", end="", file=fout)
                column_id += 1
        print(f"{column_id:13d}", end="", file=fout)
        print(file=fout)

        print(delimiter, file=fout)

    # --------------------------------------------------------------------------

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
            bead_type = 'U'
            if bead.type == "permanent":
                bead_type = 'P'
            elif bead.type == "normal":
                bead_type = 'N'
            elif bead.type == "free":
                bead_type = 'F'
            mode = "--"
            status = "--"
            client_id = "--"
            kangle = bead.kangle

            alpha = 0.0 if bead.Alpha is None else float(bead.Alpha)
            d_ad_alpha = 0.0 if bead.dAdAlpha is None else float(bead.dAdAlpha)
            free_energy_surf = 0.0 if bead.Asurf is None else float(bead.A)
            free_energy_int = 0.0 if bead.A is None else float(bead.Asurf)

            print(
                f"  {bead_id:4d} {bead_type:>6} {mode:>2} {status:>2} "
                f"{kangle:5.1f} "
                f"{alpha:6.4f} "
                f"{d_ad_alpha:12.5e} "
                f"{free_energy_surf:12.5e} "
                f"{client_id:>7}"
                f"{self.STMStep:8d}",
                end="",
                file=fout,
            )

            # CV coordinates in physical units.
            for i, cv in enumerate(self.cvs):
                unscaled = cv.unscale(bead.Pos[i])
                print(f" {float(unscaled):12.5e}", end="", file=fout)

            # CV coordinates in scaled units.
            for i, cv in enumerate(self.cvs):
                print(f" {float(bead.Pos[i]):12.5e}", end="", file=fout)

            # dA/dCV in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.Grad[i]):12.5e}", end="", file=fout)

            # dCV/dalpha in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.dCVdAlpha[i]):12.5e}", end="", file=fout)

            # Projected force, i.e. perpendicular component in scaled coordinates.
            for i in range(self.ncvs):
                print(f" {float(bead.pGrad[i]):12.5e}", end="", file=fout)

            print(f" {free_energy_int:12.5e}", end="", file=fout)

            print(file=fout)

    # --------------------------------------------------------------------------

    def write_trajectory_header(self, fout=None):
        
        if fout is None:
            fout = sys.stdout

        print(f"# STMTRAJ {self.ncvs} {self.nbeads}", file=fout)
        self.print_path_summary_header(fout)

    # --------------------------------------------------------------------------

    def write_trajectory_snapshot(self, fout=None):
        
        if fout is None:
            fout = sys.stdout

        print(f"# STMSNAP {self.STMStep}", file=fout)
        self.print_path_summary_data(fout)
        print("",file=fout) # necessary for gnuplot

# ==============================================================================
# path statistics
# ==============================================================================

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

    # --------------------------------------------------------------------------

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

    # --------------------------------------------------------------------------

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
        for bidx, bead in enumerate(self.beads,start=1):
            if bead.type == "permanent":
                continue
            
            bmov = float(np.linalg.norm(bead.Pos - bead.OPos))
            mfsize = float(np.linalg.norm(bead.uGrad))

            self.AveBeadMove += bmov
            if bmov > self.MaxBeadMove:
                self.MaxBeadMove = bmov
                self.MaxBeadMoveID = bidx

            self.AvepMFSize += mfsize
            if mfsize > self.MaxpMFSize:
                self.MaxpMFSize = mfsize
                self.MaxpMFSizeID = bidx

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
        self.MABufPLenChange[-1] = abs(self.PLenChange)
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

    # --------------------------------------------------------------------------

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

# ==============================================================================
# CLI
# ==============================================================================

def parse_args():
    
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

    print("")
    parser = argparse.ArgumentParser(
            description="The simplified string-method (STM) path optimisation on a 2D FES/PES represented by Gaussian RBFs."
        )

    # -------------------------------------------------------------------------
    # CV1
    # -------------------------------------------------------------------------

    cv1group = parser.add_argument_group("The first collective variable (CV1) specification")

    cv1group.add_argument( "--cv1label", type=str, default=r"cv1",
        help="Label of the first collective variable." )

    cv1group.add_argument( "--cv1min", type=float, required=True,
        help="Minimum value of the first collective variable." )

    cv1group.add_argument( "--cv1max", type=float, required=True,
        help="Maximum value of the first collective variable." )

    cv1group.add_argument( "--cv1nbins", type=int, required=True,
        help="Number of grid bins/points for the first collective variable." )

    # -------------------------------------------------------------------------
    # CV2
    # -------------------------------------------------------------------------

    cv2group = parser.add_argument_group("The second collective variable (CV2) specification")

    cv2group.add_argument( "--cv2label", type=str, default=r"cv2",
        help="Label of the second collective variable." )

    cv2group.add_argument( "--cv2min", type=float, required=True,
        help="Minimum value of the second collective variable." )

    cv2group.add_argument( "--cv2max", type=float, required=True,
        help="Maximum value of the second collective variable." )

    cv2group.add_argument( "--cv2nbins", type=int, required=True,
        help="Number of grid bins/points for the second collective variable." )

    # -------------------------------------------------------------------------
    # Energy label
    # -------------------------------------------------------------------------

    enegroup = parser.add_argument_group("The energy axis specification")

    enegroup.add_argument( "--enelabel", type=str, default=r"${\Delta}G [kcal/mol]$",
        help="Energy label." )

    enegroup.add_argument( "--zmax", type=float, required=True,
        help="Maximum energy value considered." )

    enegroup.add_argument( "--contour_spacing", type=float, default=1.0,
        help="Contour spacing." )

    # -------------------------------------------------------------------------
    # RBF interpolation
    # -------------------------------------------------------------------------

    rbfgroup = parser.add_argument_group("The RBF (Radial Basis Function) interpolation specification")

    rbfgroup.add_argument( "--cv1nrbfs", type=int, default=20,
        help="Number of RBFs for the first collective variable." )

    rbfgroup.add_argument( "--cv2nrbfs", type=int, default=20,
        help="Number of RBFs for the second collective variable." )

    rbfgroup.add_argument( "--rbfsx", type=float, default=1.5,
        help="Width factor for CV1 in the RBF static width mode." )

    rbfgroup.add_argument( "--rbfsy", type=float, default=1.5,
        help="Width factor for CV2 in the RBF static width mode." )

    rbfgroup.add_argument( "--rcond", type=float, default=1.0e-9,
        help="SVD cutoff for RBF fitting." )

    # -------------------------------------------------------------------------
    # Files
    # -------------------------------------------------------------------------

    filegroup = parser.add_argument_group("The input/output files specification")

    filegroup.add_argument( "--input-fes", type=str, required=True,
        help="Input FES/PES file. First three columns are CV1, CV2, energy." )
    
    filegroup.add_argument("--input-fes-x-column", type=int, default=1,
        help="Index of x-column in the input FES file." )
    
    filegroup.add_argument("--input-fes-y-column", type=int, default=2,
        help="Index of y-column in the input FES file." )
    
    filegroup.add_argument("--input-fes-e-column", type=int, default=3,
        help="Index of e-column in the input FES file." )

    filegroup.add_argument( "--input-path", type=str, required=True,
        help="Input path in the PMFLib format to optimize." )

    filegroup.add_argument( "--output-path", type=str, default="_stm.path", 
        help="Path in the PMFLib format, printed at the beginning and end of STM optimization." )
    
    filegroup.add_argument( "--summary", type=str, default="_stm.results", 
        help="Path summary, printed at the beginning and end of STM optimization." )

    filegroup.add_argument( "--trajectory", type=str, default="_stm.traj", 
        help="Path summary trajectory." )

    filegroup.add_argument( "--optlog", type=str, default="_stm.log", 
        help="STM optimization log file." )
 
    # -------------------------------------------------------------------------
    # Path
    # -------------------------------------------------------------------------
    
    pathgroup = parser.add_argument_group("Path specification")

    pathgroup.add_argument("--cvspline", type=int, default=1,
        help="Type of CV spline: 0 - interpolating cubic spline, 1 - smoothing cubic spline" )

    pathgroup.add_argument("--spline-lambda", type=float, default=0.00000002,
        help="Lambda for the internal smoothing cubic spline; 0.0 gives interpolation." )
        
    # -------------------------------------------------------------------------
    # STM Setup
    # -------------------------------------------------------------------------

    stmgroup = parser.add_argument_group("Path specification")

    stmgroup.add_argument("--nstepmax", type=int, default=500, 
        help="Maximum optimisation steps." )
    
    stmgroup.add_argument("--sfac", type=float, default=0.0,
        help="Path smoothing factor." )
    
    stmgroup.add_argument("--smoothinterval", type=int, default=0,
        help="How often to smooth the path." )
    
    stmgroup.add_argument("--reparaminterval", type=int, default=1,
        help="How often to reparametrize the path." )
    
    stmgroup.add_argument("--detect-minima-at-step", type=int, default=0,
        help="Detect minima along the pathway at given STM optimization step." )
    
    stmgroup.add_argument("--minima-threshold", type=float, default=0.5,
        help="Minimum energy separating minima along the pathway." )
    
    # -------------------------------------------------------------------------
    # Optimizer Setup
    # -------------------------------------------------------------------------

    adagroup = parser.add_argument_group("Optimizer specification")

    adagroup.add_argument("--stepsize", type=float, default=0.0002,
        help="Optimisation time step." )
    
    stmgroup.add_argument("--scale-lambda", type=float, default=10.0,
        help="Factor scaling the smoothing cubic spline lambda when minima detected." )
    
    stmgroup.add_argument("--scale-stepsize", type=float, default=0.5,
        help="Factor scaling the step size when minima detected." )

    # -------------------------------------------------------------------------
    # Termination
    # -------------------------------------------------------------------------

    termgroup = parser.add_argument_group("Termination criteria for the STM path optimization")

    termgroup.add_argument("--mabuflen", type=int, default=3,
        help="Moving-average buffer length." )

    termgroup.add_argument("--final-plenchange", type=float, default=0.001,
        help="Final threshold for moving-average path-length change." )
    
    termgroup.add_argument("--final-maxbeadmove", type=float, default=0.005,
        help="Final threshold for moving-average maximum bead movement." )

    termgroup.add_argument("--final-avebeadmove", type=float, default=.005,
        help="Final threshold for moving-average average bead movement." )
    
    termgroup.add_argument("--final-maxpmfsize", type=float, default=5.00,
        help="Final threshold for moving-average maximum projected mean-force size." )

    termgroup.add_argument("--final-avepmfsize", type=float, default=1.00,
        help="Final threshold for moving-average average projected mean-force size." )

    # -------------------------------------------------------------------------
    # Plots
    # -------------------------------------------------------------------------

    plotgroup = parser.add_argument_group("The graphical plot specification")

    plotgroup.add_argument("--plot", action="store_true", default=False,
        help="Write PNG plots." )
    
    plotgroup.add_argument("--plot-prefix", type=str, default="_stm", 
        help="Prefix for output PNG plots." )
    
    plotgroup.add_argument("--show", action="store_true", default=False,
        help="Show plots interactively after saving." )

    plotgroup.add_argument('--figsize',type=parse_figsize,default=(6.4, 4.8),  # Default Matplotlib size fallback
        help="Figure size as 'width,height' in inches (default: 6.4,4.8)" )

    plotgroup.add_argument( "--dpi", type=int, default=300,
        help="Resolution for plot figures." )

    # -------------------------------------------------------------------------

    args = parser.parse_args()

    # check some input 
    if args.mabuflen < 1:
        parser.error("--mabuflen must be at least 1")

    if not (0.0 <= args.sfac <= 1.0):
        parser.error("--sfac must be in the interval [0, 1]")

    if args.cvspline not in (0, 1):
        parser.error("--cvspline must be 0 or 1")

    if args.minima_threshold < 0.0:
        parser.error("--minima-threshold > 0.")

    return args

# ==============================================================================
# Main
# ==============================================================================

def main() -> None:

    print("")
    print("#==============================================================================#")
    print("#          *** Simplified String Method (STM) on 2D Energy Surface ***         #")
    print("#         The stm-path-2D-surface utility is part of the PMFLib toolkit.       #")
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

    print("# All arguments ...")
    print(options)

    # do all STM stuff :-)
    stmpath = STMPath(args)
    stmpath.stm_optimize(args)

    print("")
    
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

# ==============================================================================
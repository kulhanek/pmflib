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

# ==============================================================================

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

# ==============================================================================

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
        # BUG??? self.sd[1] = self.y[1] - mu*(f[1]*q[1] + r[1]*q[2])*self.sigma[0]
        self.sd[1] = self.y[1] - mu*(f[1]*q[1] + r[1]*q[2])*self.sigma[1]
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


# ==============================================================================

class CVSmoothingCubicSplineSVD:
    """
    Natural smoothing cubic spline solved by SVD.

    The spline minimizes

        sum_i w_i * (y_i - z_i)^2
        + lam * integral (S''(x))^2 dx

    where

        w_i = 1 / sigma_i^2

    Interface
    ---------
    spline = CVSmoothingCubicSplineSVD(lam, all_sigma, rcond)

    spline.update_points(x, y)

    value = spline(x)
    der1  = spline(x, 1)
    der2  = spline(x, 2)
    der3  = spline(x, 3)

    Notes
    -----
    lam = 0.0 gives a natural interpolating cubic spline.
    Larger lam gives stronger smoothing.
    """

    def __init__(self, lam: float, all_sigma: float | np.ndarray, rcond: float):
        self.lam = float(lam)
        self.all_sigma = all_sigma
        self.rcond = float(rcond)

        if self.lam < 0.0:
            raise ValueError("lam must be non-negative.")

        if self.rcond < 0.0:
            raise ValueError("rcond must be non-negative.")

        self.x = None
        self.y = None
        self.sigma = None
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

        sigma = self._prepare_sigma(x.size)
        sigma = sigma[order]

        if np.any(sigma <= 0.0):
            raise ValueError("All sigma values must be positive.")

        h = np.diff(x)

        if np.any(h <= 0.0):
            raise ValueError("x values must be strictly increasing.")

        self.x = x
        self.y = y
        self.sigma = sigma
        self.w = 1.0 / sigma**2

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
            raise RuntimeError("Spline is not initialized. Call update_points(x, y) first.")

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
        """

        q, r = self._build_reinsch_matrices()

        # Compute K = Q R^{-1} Q^T.
        #
        # We do not explicitly invert R. Instead, solve
        #
        #     R X = Q^T
        #
        # by SVD.
        rinv_qt = self._svd_solve(r, q.T)
        k = q @ rinv_qt

        lhs = np.diag(self.w) + self.lam * k
        rhs = self.w * self.y

        z = self._svd_solve(lhs, rhs)

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

        m_inner = self._svd_solve(amat, rhs)

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
    # Helpers
    # -------------------------------------------------------------------------

    def _prepare_sigma(self, npts):
        """
        Prepare sigma array from self.all_sigma.

        self.all_sigma can be either a scalar or an array with length npts.
        """

        sigma = np.asarray(self.all_sigma, dtype=float)

        if sigma.ndim == 0:
            return np.full(npts, float(sigma), dtype=float)

        if sigma.ndim != 1:
            raise ValueError("all_sigma must be either a scalar or a one-dimensional array.")

        if sigma.size != npts:
            raise ValueError("If all_sigma is an array, it must have the same length as x and y.")

        return sigma.copy()

    def _svd_solve(self, a, b):
        """
        Solve a*x = b using SVD pseudoinverse.

        Singular values are accepted if

            s_i > rcond * max(s)

        Parameters
        ----------
        a : ndarray
            Matrix.

        b : ndarray
            Right-hand side vector or matrix.
        """

        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)

        u, s, vt = np.linalg.svd(a, full_matrices=False)

        if s.size == 0:
            raise np.linalg.LinAlgError("SVD failed: no singular values found.")

        cutoff = self.rcond * np.max(s)

        sinv = np.zeros_like(s)
        keep = s > cutoff
        sinv[keep] = 1.0 / s[keep]

        if b.ndim == 1:
            return vt.T @ (sinv * (u.T @ b))

        return vt.T @ (sinv[:, None] * (u.T @ b))


# class CVSmoothingCubicSplineSVD:
#     """
#     Natural smoothing cubic spline using SVD-based linear solves.

#     The spline minimizes

#         lam * sum_i ((y_i - z_i)^2 / sigma_i^2)
#         + (1.0 - lam) * integral (S''(x))^2 dx

#     Interface
#     ---------
#     spline = CVSmoothingCubicSplineSVD(lam, all_sigma, rcond)

#     spline.update_points(x, y)

#     value = spline(x)
#     der1  = spline(x, 1)
#     der2  = spline(x, 2)
#     der3  = spline(x, 3)

#     Parameters
#     ----------
#     lam : float
#         Smoothing/data-balance parameter in [0, 1].

#         lam = 1:
#             Interpolating natural cubic spline through y.

#         lam = 0:
#             Pure roughness minimization. The minimizer is not unique because
#             any straight line has zero roughness. With the SVD pseudoinverse,
#             the minimum-norm solution is selected.

#     all_sigma : float or array_like
#         If scalar, the same sigma is used for all points.
#         If array-like, it must have the same length as x and y.

#     rcond : float
#         Relative singular-value cutoff for SVD pseudoinverse solves.
#     """

#     def __init__(self, lam, all_sigma, rcond):
#         self.lam = float(lam)
#         self.all_sigma = all_sigma
#         self.rcond = float(rcond)

#         if self.lam < 0.0 or self.lam > 1.0:
#             raise ValueError("lam must be in the interval [0, 1].")

#         if self.rcond < 0.0:
#             raise ValueError("rcond must be non-negative.")

#         self.x = None
#         self.y = None
#         self.sigma = None

#         self.n = 0
#         self.h = None

#         self.z = None

#         self.a = None
#         self.b = None
#         self.c = None
#         self.d = None

#     # -------------------------------------------------------------------------
#     # Public interface
#     # -------------------------------------------------------------------------

#     def merge_close_points(self, x, y, sigma=None, atol=1e-12):
#         """
#         Merge duplicate or nearly duplicate x values.

#         Points with distance <= atol are treated as one point.
#         Their x, y, and optionally sigma values are averaged.

#         Parameters
#         ----------
#         x : array_like
#             Input x values.

#         y : array_like
#             Input y values.

#         sigma : array_like or None
#             Optional sigma values.

#         atol : float
#             Absolute tolerance for considering two x values identical.

#         Returns
#         -------
#         x_new, y_new, sigma_new
#             Merged arrays. sigma_new is None if sigma was None.
#         """

#         x = np.asarray(x, dtype=float)
#         y = np.asarray(y, dtype=float)

#         if sigma is not None:
#             sigma = np.asarray(sigma, dtype=float)

#         order = np.argsort(x)
#         x = x[order]
#         y = y[order]

#         if sigma is not None:
#             sigma = sigma[order]

#         x_new = []
#         y_new = []
#         sigma_new = [] if sigma is not None else None

#         group_x = [x[0]]
#         group_y = [y[0]]
#         group_sigma = [sigma[0]] if sigma is not None else None

#         for i in range(1, x.size):
#             if abs(x[i] - group_x[-1]) <= atol:
#                 group_x.append(x[i])
#                 group_y.append(y[i])
#                 if sigma is not None:
#                     group_sigma.append(sigma[i])
#             else:
#                 x_new.append(np.mean(group_x))
#                 y_new.append(np.mean(group_y))

#                 if sigma is not None:
#                     sigma_new.append(np.mean(group_sigma))

#                 group_x = [x[i]]
#                 group_y = [y[i]]
#                 group_sigma = [sigma[i]] if sigma is not None else None

#         x_new.append(np.mean(group_x))
#         y_new.append(np.mean(group_y))

#         if sigma is not None:
#             sigma_new.append(np.mean(group_sigma))

#         x_new = np.asarray(x_new, dtype=float)
#         y_new = np.asarray(y_new, dtype=float)

#         if sigma is not None:
#             sigma_new = np.asarray(sigma_new, dtype=float)

#         return x_new, y_new, sigma_new

#     def update_points(self, x, y):
#         x = np.asarray(x, dtype=float)
#         y = np.asarray(y, dtype=float)

#         if x.ndim != 1:
#             raise ValueError("x must be a one-dimensional array.")

#         if y.ndim != 1:
#             raise ValueError("y must be a one-dimensional array.")

#         if x.size != y.size:
#             raise ValueError("x and y must have the same length.")

#         if x.size < 2:
#             raise ValueError("At least two points are required.")

#         sigma = self._prepare_sigma(x.size)

#         x, y, sigma = self.merge_close_points(
#             x,
#             y,
#             sigma=sigma,
#             atol=1e-12,
#         )

#         if x.size < 2:
#             raise ValueError("At least two distinct x points are required.")

#         h = np.diff(x)

#         if np.any(h <= 0.0):
#             raise ValueError("x values must be strictly increasing after merging.")

#         if np.any(sigma <= 0.0):
#             raise ValueError("All sigma values must be positive.")

#         self.x = x
#         self.y = y
#         self.sigma = sigma

#         self.n = x.size
#         self.h = h

#         self._build_spline()

    # def update_points(self, x, y):
    #     """
    #     Set or update spline points.

    #     Parameters
    #     ----------
    #     x : array_like
    #         Knot positions.
    #     y : array_like
    #         Values at knot positions.
    #     """

    #     x = np.asarray(x, dtype=float)
    #     y = np.asarray(y, dtype=float)

    #     if x.ndim != 1:
    #         raise ValueError("x must be a one-dimensional array.")

    #     if y.ndim != 1:
    #         raise ValueError("y must be a one-dimensional array.")

    #     if x.size != y.size:
    #         raise ValueError("x and y must have the same length.")

    #     if x.size < 2:
    #         raise ValueError("At least two points are required.")

    #     order = np.argsort(x)
    #     x = x[order]
    #     y = y[order]

    #     h = np.diff(x)

    #     if np.any(h <= 0.0):
    #         raise ValueError("x values must be strictly increasing.")

    #     sigma = self._prepare_sigma(x.size)

    #     if sigma.size == x.size:
    #         sigma = sigma[order]

    #     if np.any(sigma <= 0.0):
    #         raise ValueError("All sigma values must be positive.")

    #     self.x = x
    #     self.y = y
    #     self.sigma = sigma

    #     self.n = x.size
    #     self.h = h

    #     self._build_spline()

    # def __call__(self, x_eval, der=0):
    #     """
    #     Evaluate spline or its derivative.

    #     Parameters
    #     ----------
    #     x_eval : float or array_like
    #         Evaluation point or points.

    #     der : int, default=0
    #         Derivative order.

    #         der = 0:
    #             value

    #         der = 1:
    #             first derivative

    #         der = 2:
    #             second derivative

    #         der = 3:
    #             third derivative

    #         der > 3:
    #             zero

    #     Returns
    #     -------
    #     float or ndarray
    #         Spline value or derivative.
    #     """

    #     if self.x is None:
    #         raise RuntimeError("Spline points are not initialized. Call update_points(x, y) first.")

    #     if der < 0:
    #         raise ValueError("Derivative order must be non-negative.")

    #     scalar_input = np.isscalar(x_eval)
    #     x_eval = np.asarray(x_eval, dtype=float)

    #     idx = np.searchsorted(self.x, x_eval, side="right") - 1
    #     idx = np.clip(idx, 0, self.n - 2)

    #     t = x_eval - self.x[idx]

    #     a = self.a[idx]
    #     b = self.b[idx]
    #     c = self.c[idx]
    #     d = self.d[idx]

    #     if der == 0:
    #         out = a + b*t + c*t**2 + d*t**3
    #     elif der == 1:
    #         out = b + 2.0*c*t + 3.0*d*t**2
    #     elif der == 2:
    #         out = 2.0*c + 6.0*d*t
    #     elif der == 3:
    #         out = 6.0*d
    #     else:
    #         out = np.zeros_like(x_eval, dtype=float)

    #     if scalar_input:
    #         return float(out)

    #     return out

    # # -------------------------------------------------------------------------
    # # Spline construction
    # # -------------------------------------------------------------------------

    # def _build_spline(self):
    #     """
    #     Build smoothing spline.
    #     """

    #     if self.n == 2:
    #         self.z = self.y.copy()
    #         self._build_natural_cubic(self.z)
    #         return

    #     if self.lam == 1.0:
    #         self.z = self.y.copy()
    #     else:
    #         self.z = self._calculate_smoothed_values()

    #     self._build_natural_cubic(self.z)

    # def _calculate_smoothed_values(self):
    #     """
    #     Calculate smoothed knot values z.
    #     """

    #     q, r = self._build_qr_matrices()

    #     rinv_qt = self._svd_solve(r, q.T)

    #     k = q @ rinv_qt

    #     w_diag = 1.0 / (self.sigma**2)

    #     lhs = self.lam * np.diag(w_diag) + (1.0 - self.lam) * k
    #     rhs = self.lam * w_diag * self.y

    #     z = self._svd_solve(lhs, rhs)

    #     return z

    # def _build_qr_matrices(self):
    #     """
    #     Build Reinsch Q and R matrices for a natural cubic spline.
    #     """

    #     n = self.n
    #     h = self.h

    #     q = np.zeros((n, n - 2), dtype=float)

    #     for i in range(n - 2):
    #         q[i, i] = 1.0 / h[i]
    #         q[i + 1, i] = -1.0 / h[i] - 1.0 / h[i + 1]
    #         q[i + 2, i] = 1.0 / h[i + 1]

    #     r = np.zeros((n - 2, n - 2), dtype=float)

    #     for i in range(n - 2):
    #         r[i, i] = (h[i] + h[i + 1]) / 3.0

    #     for i in range(n - 3):
    #         r[i, i + 1] = h[i + 1] / 6.0
    #         r[i + 1, i] = h[i + 1] / 6.0

    #     return q, r

    # def _build_natural_cubic(self, z):
    #     """
    #     Build natural cubic interpolation coefficients through smoothed values z.

    #     On interval [x_i, x_{i+1}], the spline is

    #         S_i(x) = a_i + b_i t + c_i t^2 + d_i t^3

    #     with

    #         t = x - x_i
    #     """

    #     n = self.n
    #     h = self.h

    #     if n == 2:
    #         self.a = np.array([z[0]], dtype=float)
    #         self.b = np.array([(z[1] - z[0]) / h[0]], dtype=float)
    #         self.c = np.array([0.0], dtype=float)
    #         self.d = np.array([0.0], dtype=float)
    #         return

    #     amat = np.zeros((n - 2, n - 2), dtype=float)
    #     rhs = np.zeros(n - 2, dtype=float)

    #     for i in range(1, n - 1):
    #         row = i - 1

    #         if i > 1:
    #             amat[row, row - 1] = h[i - 1]

    #         amat[row, row] = 2.0 * (h[i - 1] + h[i])

    #         if i < n - 2:
    #             amat[row, row + 1] = h[i]

    #         rhs[row] = 6.0 * (
    #             (z[i + 1] - z[i]) / h[i]
    #             - (z[i] - z[i - 1]) / h[i - 1]
    #         )

    #     m_inner = self._svd_solve(amat, rhs)

    #     m = np.zeros(n, dtype=float)
    #     m[1:-1] = m_inner

    #     self.a = z[:-1].copy()

    #     self.b = (
    #         (z[1:] - z[:-1]) / h
    #         - h * (2.0*m[:-1] + m[1:]) / 6.0
    #     )

    #     self.c = m[:-1] / 2.0

    #     self.d = (m[1:] - m[:-1]) / (6.0*h)

    # # -------------------------------------------------------------------------
    # # Helpers
    # # -------------------------------------------------------------------------

    # def _prepare_sigma(self, n):
    #     """
    #     Prepare sigma array.
    #     """

    #     sigma = np.asarray(self.all_sigma, dtype=float)

    #     if sigma.ndim == 0:
    #         return np.full(n, float(sigma), dtype=float)

    #     if sigma.ndim != 1:
    #         raise ValueError("all_sigma must be either a scalar or a one-dimensional array.")

    #     if sigma.size != n:
    #         raise ValueError("If all_sigma is an array, it must have the same length as x and y.")

    #     return sigma.copy()

    # def _svd_solve(self, a, b):
    #     """
    #     Solve a x = b by SVD pseudoinverse.

    #     Singular values are accepted if

    #         s_i > rcond * max(s)

    #     Parameters
    #     ----------
    #     a : ndarray
    #         Matrix.

    #     b : ndarray
    #         Right-hand side.

    #     Returns
    #     -------
    #     ndarray
    #         Pseudoinverse solution.
    #     """

    #     a = np.asarray(a, dtype=float)
    #     b = np.asarray(b, dtype=float)

    #     u, s, vt = np.linalg.svd(a, full_matrices=False)

    #     if s.size == 0:
    #         raise np.linalg.LinAlgError("SVD failed: no singular values found.")

    #     cutoff = self.rcond * np.max(s)

    #     sinv = np.zeros_like(s)
    #     keep = s > cutoff
    #     sinv[keep] = 1.0 / s[keep]

    #     if b.ndim == 1:
    #         return vt.T @ (sinv * (u.T @ b))

    #     return vt.T @ (sinv[:, None] * (u.T @ b))

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

        dx = 1.0 / self.x_axis.npts
        dy = 1.0 / self.y_axis.npts

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
                            path_types: np.ndarray = None,
                            filename: str = None, show: bool = None,
                            figsize = None, dpi: int = 300, legend: bool = False) -> None:

            cmap, norm = self.make_colormap()
            cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing
            levels = np.arange(0.0, cmax + self.contour_spacing, self.contour_spacing)

            bead_colors = {
                "flexible":  "green",
                "terminal":  "green",
                "permanent": "black",
                "kink":      "orange",
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
                for bead_type in ("flexible", "terminal", "permanent", "kink"):
                    mask = path_types == bead_type
                    if not np.any(mask):
                        continue
                    ax.scatter(
                        path_x[mask], path_y[mask],
                        s=28,
                        c=bead_colors[bead_type],
                        edgecolors="white",
                        linewidths=0.6,
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

        self.type       = None

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

    def UpdatePositionAdaBelief(self,step,beta1,beta2,mingnormeps,cvs):

        if self.type == "permanent":
            return
        
        grad = np.zeros(self.ncvs)

        if self.type == "terminal" or self.type == "kink":
            grad[:] = self.Grad[:]  # gradient descent move
        else:
            grad[:] = self.pGrad[:] # perpedicular move

        self.mt[:] = beta1 * self.mt[:] + (1.0 - beta1) * grad[:]
        self.vt[:] = beta2 * self.vt[:] + (1.0 - beta2) * ((grad[:] - self.mt[:])*(grad[:] - self.mt[:]) + mingnormeps)

        self.beta1t = self.beta1t * beta1
        self.beta2t = self.beta2t * beta2

        for i in range(self.ncvs):

            mthat = self.mt[i]/(1.0-self.beta1t)
            vthat = self.vt[i]/(1.0-self.beta2t)

            dm = step * mthat/(math.sqrt(vthat)+mingnormeps)

            if (cvs[i].smaxmov <= 0) or (math.fabs(dm) < cvs[i].smaxmov):
                self.Pos[i] = self.Pos[i] - dm
            else:
                self.Pos[i] = self.Pos[i] - cvs[i].smaxmov*math.copysign(1.0,dm)

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
        self.surface.load(args.input_fes)

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
            print(f"      Lambda: {args.spline_lambda:10.5f}")
            print(f"      Sigma:  {args.spline_sigma:10.5f}")
            self.cv_splines =  [
                CVSplineSmoothingCubic(lam=args.spline_lambda, sigma=args.spline_sigma)
                for _ in range(self.ncvs)
            ]
        elif args.cvspline == 2:
            print(f"  >>> Smoothing Cubic Spline via SVD")
            print(f"      Lambda: {args.spline_lambda:10.5f}")
            print(f"      RCond:  {args.spline_rcond:10.6e}")
            self.cv_splines =  [
                CVSmoothingCubicSplineSVD(lam=args.spline_lambda, all_sigma=1.0, rcond=args.spline_rcond)
                for _ in range(self.ncvs)
            ]
        else:
            print(f"  >>> Interpolating Cubic Spline")
            self.cv_splines =  [CVSplineInterpolatingCubic() for _ in range(self.ncvs)]

        self.PathParamMode = args.path_param_mode

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

        # setup AdaBelief
        self.StepSize           = args.stepsize
        self.AdamB1             = args.beta1
        self.AdamB2             = args.beta2
        self.MinGNormEps        = args.mingnormesp

        # smoothing, reparameterization
        self.SmoothInterval     = args.smoothinterval
        self.SmoothingFac       = args.sfac
        self.ReparamInterval    = args.reparaminterval

        # STM
        self.STMStep            = 0
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

        self.PathParamMode :  0 - "chord-length" parameterization
                              1 - centripetal parameterization

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
        total_length_sqrt = 0.0

        for b in range(1, len(beads)):
            diff = beads[b].Pos - beads[b - 1].Pos
            slen = np.linalg.norm(diff)

            total_length      += slen
            total_length_sqrt += math.sqrt(slen)

        if total_length == 0.0:
            raise RuntimeError("path has zero length")

        # ---------------------------------------------------------------------
        # Initial alpha values from linear interpolation.
        # ---------------------------------------------------------------------

        beads[0].Alpha = 0.0

        path_length = 0.0
        path_length_sqrt = 0.0

        for b in range(1, len(beads) - 1):
            diff = beads[b].Pos - beads[b - 1].Pos
            slen = np.linalg.norm(diff)

            if slen == 0.0:
                raise RuntimeError("path segment has zero length")
            
            path_length      += slen
            path_length_sqrt += math.sqrt(slen)

            if self.PathParamMode == 1:
                beads[b].Alpha = path_length_sqrt / total_length_sqrt
            else:
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

    def update_all_positions(self):

        # backup old positions
        for bead in self.beads:
            bead.OPos[:] = bead.Pos[:]

        self.OldCPathLength = self.CurrentPathLength

        self.STMStep += 1

        # update positions
        for bead in self.beads:
            # bead types are handled in UpdatePositionAdaBelief()
            bead.UpdatePositionAdaBelief(self.StepSize,self.AdamB1,self.AdamB2,self.MinGNormEps,self.cvs)
            
    # --------------------------------------------------------------------------

    def smooth_all_positions(self):

        if (self.SmoothInterval == 0) or (self.STMStep % self.SmoothInterval != 0 ):
            return;
    
        old_pos = np.array([bead.Pos.copy() for bead in self.beads])

        for i in range(1, self.nbeads - 1):
            if self.beads[i].type == "permanent" or self.beads[i].type == "kink":
                # skip kink or permanent beads 
                continue

            self.beads[i].Pos[:] = (
                (1.0 - self.SmoothingFac) * old_pos[i]
                + 0.5 * self.cSmoothingFac * (old_pos[i - 1] + old_pos[i + 1])
            )

    # --------------------------------------------------------------------------

    def reparametrize_all_positions(self):

        if (self.ReparamInterval == 0) or (self.STMStep % self.ReparamInterval != 0 ):
            return;
    
        self.parametrize_path(self.beads)

        # generate evenly distributed set of alphas
        for i in range(self.nbeads):
            self.beads[i].Alpha = float(i) / float(self.nbeads-1)

        self.beads[0].Alpha = 0.0
        self.beads[-1].Alpha = 1.0

        # update positions along the path based on new alphas
        for cv in range(self.ncvs):
            # FIXME
            # redistribute other beads
            for bead in self.beads:
                if bead.type == "permanent":
                    continue
                bead.Pos[cv] = self.cv_splines[cv](bead.Alpha)
       
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

        self.CurrentPathLength = self.parametrize_path(self.beads)

        for bead in self.beads:            
            # Calculate MF and A
            bead.Asurf, grad, hessian = self.surface.eval_uv(bead.Pos)
            bead.Grad[:] = grad[:]

            # Calculate tangent dCV/dalpha from the path splines.
            for i in range(self.ncvs):
                bead.dCVdAlpha[i] = self.cv_splines[i](bead.Alpha, 1)

            slen2 = float(np.dot(bead.dCVdAlpha, bead.dCVdAlpha))
            if slen2 == 0.0:
                raise RuntimeError("derivative segment has zero length")

            # Projector perpendicular to the path:
            #     P = I - t t^T / |t|^2

            bead.P[:, :] = np.eye(self.ncvs) - np.outer(bead.dCVdAlpha, bead.dCVdAlpha) / slen2

        # Project gradients
        for bead in self.beads:
            bead.pGrad[:] = bead.P @ bead.Grad

        v1 = np.zeros(self.ncvs)
        v2 = np.zeros(self.ncvs)

        # Calculate kink angles
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
            raise RuntimeError("At least two input beads (flexible/permanent/terminal/kink) must be specified in the input PATH file!")
        
        if self.nbeads < 2:
            raise RuntimeError("At least two beads (nbeads) must be requested in the input PATH file!")

        if input_beads[0].type == "kink":
            raise RuntimeError(f"Terminal bead must be terminal/flexible/permanent but {input_beads[0].type} was specified!")
        if input_beads[-1].type == "kink":
            raise RuntimeError(f"Terminal bead must be terminal/flexible/permanent but {input_beads[0].type} was specified!")

        if len(input_beads) == self.nbeads:
            # check boundaries
            self.check_boundaries_of_beads(input_beads)

            # path is complete, keep it as it is
            return input_beads
        
        # path is incompleted, it will be rebuilded

        for bead in input_beads[1:-1]:
            if bead.type != "flexible":
                raise RuntimeError(f"For the incomplete path, all inner beads must be flexible but {bead.type} was requested!")

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
            bead.type  = "flexible"

            # copy type of terminals from the input path
            if (bidx == 0) or (bidx == self.nbeads - 1):
                bead.type = input_beads[0].type

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
        permanent
        terminal
        flexible
        kink
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

                    elif key in ("permanent", "terminal", "flexible", "kink"):
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
            raise ValueError(f"at least two beads must be provided | flexible/permanent/terminal/kink")
        
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
                # Python version supports only flexible beads.
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
        print("#  ID   Type  MO ST KinkA  alpha    dA/dalpha            A     CID Updates", end="", file=fout)
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

        print(f"          ENE", end="", file=fout)
        print(file=fout)

        # Delimiters.
        delimiter = "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------"
        delimiter += " ------------" * (5 * self.ncvs)
        delimiter += " ------------"
        print(delimiter, file=fout)

        # CV metadata rows.
        print(f"{'#      names':<68}", end="", file=fout)
        for _ in range(5):
            for cv in self.cvs:
                print(f" {cv.name:>12}", end="", file=fout)
        print(file=fout)

        print(f"{'#      types':<68}", end="", file=fout)
        for _ in range(2):
            for cv in self.cvs:
                print(f" {cv.type:>12}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      min':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.pathmin:12.5e}", end="", file=fout)
        for cv in self.cvs:
            print(f" {0.0:12.5e}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      max':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.pathmax:12.5e}", end="", file=fout)
        for cv in self.cvs:
            print(f" {1.0:12.5e}", end="", file=fout)
        for _ in range(3):
            for cv in self.cvs:
                print(f"             ", end="", file=fout)
        print(file=fout)

        print(f"{'#      maxmov':<68}", end="", file=fout)
        for cv in self.cvs:
            print(f" {cv.maxmov:12.5e}", end="", file=fout)
        for _ in range(3):
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
            elif bead.type == "flexible":
                bead_type = 'F'
            elif bead.type == "terminal":
                bead_type = 'T'
            elif bead.type == "kink":
                bead_type = 'K'
            mode = "--"
            status = "--"
            client_id = "--"
            kangle = bead.kangle

            alpha = 0.0 if bead.Alpha is None else float(bead.Alpha)
            d_ad_alpha = 0.0 if bead.dAdAlpha is None else float(bead.dAdAlpha)
            free_energy = 0.0 if bead.A is None else float(bead.A)
            ene_surf = 0.0 if bead.Asurf is None else float(bead.Asurf)

            print(
                f"  {bead_id:4d} {bead_type:>6} {mode:>2} {status:>2} "
                f"{kangle:5.1f} "
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

            print(f" {ene_surf:12.5e}", end="", file=fout)

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
            
            # print(bead.Pos,bead.OPos)
            bmov = float(np.linalg.norm(bead.Pos - bead.OPos))

            if bead.type == "flexible":
                mfsize = float(np.linalg.norm(bead.pGrad))
            else:
                # kink, terminal
                mfsize = float(np.linalg.norm(bead.Grad))

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

    enegroup.add_argument(
        "--enelabel", type=str, default=r"${\Delta}G [kcal/mol]$",
        help="Energy label."
    )

    enegroup.add_argument(
        "--zmax", type=float, required=True,
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

    pathgroup.add_argument("--cvspline", type=int, default=2,
        help="Type of CV spline: 0 - interpolating cubic spline, 1 - smoothing cubic spline, 2 - smoothing cubic spline SVD." )

    pathgroup.add_argument("--spline-lambda", type=float, default=0.999,
        help="Lambda for the internal smoothing cubic spline; 1.0 gives interpolation." )
    
    pathgroup.add_argument("--spline-sigma", type=float, default=0.01,
        help="Default sigma assigned to all knots of the internal smoothing cubic spline." )
    
    pathgroup.add_argument("--spline-rcond", type=float, default=1e-9,
        help="SVD cutoff for CVSmoothingCubicSplineSVD." )

    pathgroup.add_argument("--path-param-mode", type=int, default=1,
        help="Path parameterization mode: 0 - 'chord-length' parameterization, 1 - centripetal parameterization." )
    
    # -------------------------------------------------------------------------
    # STM Setup
    # -------------------------------------------------------------------------

    stmgroup = parser.add_argument_group("Path specification")

    stmgroup.add_argument("--nstepmax", type=int, default=200, 
        help="Maximum optimisation steps." )
    
    stmgroup.add_argument("--sfac", type=float, default=0.0,
        help="Path smoothing factor."
    )
    
    stmgroup.add_argument("--smoothinterval", type=int, default=0,
        help="How often to smooth the path." )
    
    stmgroup.add_argument("--reparaminterval", type=int, default=1,
        help="How often to reparametrize the path." )

    # -------------------------------------------------------------------------
    # STM Setup
    # -------------------------------------------------------------------------

    adagroup = parser.add_argument_group("AdaBelief specification")

    adagroup.add_argument("--stepsize", type=float, default=0.003,
        help="Optimisation time step." )

    adagroup.add_argument("--beta1", type=float, default=0.7,
        help="AdaBelief beta1." )

    adagroup.add_argument("--beta2", type=float, default=0.99,
        help="AdaBelief beta2." )

    adagroup.add_argument("--mingnormesp", type=float, default=1e-7,
        help="AdaBelief epsilon." )

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
    
    termgroup.add_argument("--final-maxpmfsize", type=float, default=2.00,
        help="Final threshold for moving-average maximum projected mean-force size." )

    termgroup.add_argument("--final-avepmfsize", type=float, default=0.80,
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

    if args.cvspline not in (0, 1, 2):
        parser.error("--cvspline must be 0 or 1 or 2")

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
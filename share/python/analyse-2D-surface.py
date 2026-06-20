#!/usr/bin/env python3

import argparse
import time
import math
import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, ListedColormap
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from skimage.segmentation import watershed

# ------------------------------------------------------------------------------

rad2deg = 180.0 / np.pi
rfac = 0.001987204258640       # kcal/mol/K

# ==============================================================================
# CV/Axis
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

        self.cvmin = float(cvmin)
        self.cvmax = float(cvmax)
        self.nrbfs = int(nrbfs)
        self.npts = int(npts)
        self.label = label

        self.range = self.cvmax - self.cvmin
        self.width = 1.0 / self.nrbfs
        self.width_scale = 1.0

        # Non-periodic CV: centers at both boundaries, therefore nx+1 centers
        self.centers = np.arange(self.nrbfs + 1, dtype=float) / self.nrbfs

    def scale(self, x):
        """Convert a physical CV value to the internal scaled coordinate."""
        u = (x - self.cvmin) / self.range
        return u

    def unscale(self, u):
        """Convert an internal scaled coordinate to the physical CV value."""
        return self.cvmin + u * self.range

    def delta(self, u):
        """
        Difference between scaled coordinate u and all RBF centers.
        """
        du = u - self.centers
        return du

# ==============================================================================
# EnergySurface
# ==============================================================================

class EnergySurface2D:
    """
    2D energy surface represented by Gaussian radial basis functions.

    Internal coordinates are scaled coordinates:

        u = (x - xmin) / (xmax - xmin)
        v = (y - ymin) / (ymax - ymin)

    The RBF expansion is

        E(u, v) = sum_i sum_j A_ij phi_ij(u, v)

    where

        phi_ij = exp[-1/2 ((u - ui)^2 / su^2 + (v - vj)^2 / sv^2)]

    The method eval(x) returns derivatives with respect to original coordinates.
    The method eval_uv(u) returns derivatives with respect to scaled coordinates.
    """

    def __init__(
        self,
        cv1min,
        cv1max,
        cv1_nx,
        cv1_px,
        cv1_label,
        cv2min,
        cv2max,
        cv2_nx,
        cv2_px,
        cv2_label,
        ene_label,
        zmax,
        contour_spacing,
        thrfac,
        temp,
        random_seed
    ):
        """
        npx, npy : int
            Number of grid points along the first and second CV.

        zmax : float
            maximum Z value

        contour_spacing : float
            Spacing between contour lines in energy units.
        """

        self.x_axis = Axis(cv1min, cv1max, cv1_nx, cv1_px, cv1_label)
        self.y_axis = Axis(cv2min, cv2max, cv2_nx, cv2_px, cv2_label)

        self.x_data = None
        self.y_data = None
        self.e_data = None


        self.X = None
        self.Y = None
        self.Z = None

        self.U = None
        self.V = None
        self.SNG = None

        self.temp = float(temp)
        self.thr = float(thrfac) * temp * rfac

        self.ene_label = ene_label
        self.zmax = float(zmax)
        self.contour_spacing = float(contour_spacing)

        self.amplitudes = None

        self.sp_guesses = []
        self.sp_optimized = []

        if random_seed is None:
            random_seed = int(time.time())
        print(f"  Random seed: {random_seed}")
        self.rng = np.random.default_rng(random_seed)

    # --------------------------------------------------------------------------
    # Coordinate conversion
    # --------------------------------------------------------------------------

    def to_scaled(self, x):
        """
        Convert original coordinates to scaled coordinates.

        Parameters
        ----------
        x : array-like with two elements
            Original coordinates [x, y].

        Returns
        -------
        u : ndarray, shape (2,)
            Scaled coordinates [u, v].
        """
        x = np.asarray(x, dtype=float)

        if x.shape != (2,):
            raise ValueError("x must be a vector with two elements: [x, y]")

        return np.array(
            [
                self.x_axis.scale(x[0]),
                self.y_axis.scale(x[1]),
            ],
            dtype=float,
        )

    # --------------------------------------------------------------------------

    def from_scaled(self, u):
        """
        Convert scaled coordinates to original coordinates.

        Parameters
        ----------
        u : array-like with two elements
            Scaled coordinates [u, v].

        Returns
        -------
        x : ndarray, shape (2,)
            Original coordinates [x, y].
        """
        u = np.asarray(u, dtype=float)

        if u.shape != (2,):
            raise ValueError("u must be a vector with two elements: [u, v]")

        return np.array(
            [
                self.x_axis.unscale(u[0]),
                self.y_axis.unscale(u[1]),
            ],
            dtype=float,
        )

    # --------------------------------------------------------------------------
    # Data loading and fitting
    # --------------------------------------------------------------------------

    def load(self, filename, xcolumn=1, ycolumn=2, ecolumn=3):
        """
        Load input data from a text file.

        Comments and empty lines are ignored.
        The first three columns are interpreted as:

            column 1: cv1
            column 2: cv2
            column 3: energy
        """
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

    # --------------------------------------------------------------------------

    def calc_ene_and_sng(self):

        if self.x_axis.npts <= 1 or self.y_axis.npts <= 1:
            raise ValueError("npx and npy must be larger than 1")

        # Grid in original coordinates
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
        self.x_grid = 0.5 * (x_edges[:-1] + x_edges[1:])
        self.y_grid = 0.5 * (y_edges[:-1] + y_edges[1:])

        self.X, self.Y = np.meshgrid(self.x_grid, self.y_grid, indexing="xy")
        self.Z = np.empty_like(self.X, dtype=float)

        # Evaluate the surface through eval(), as requested
        for iy in range(self.y_axis.npts):
            for ix in range(self.x_axis.npts):
                value, gradient, hessian = self.eval([self.X[iy, ix], self.Y[iy, ix]])
                self.Z[iy, ix] = value

        # detect unsampled points
        self.unsampled_mask, nearest_dist = self.detect_unsampled_grid_points()

        # set them to zmax
        self.Z[self.unsampled_mask] = self.zmax

        # Grid in scaled coordinates
        u_edges = np.linspace(
            0.0,
            1.0,
            self.x_axis.npts + 1,
        )
        v_edges = np.linspace(
            0.0,
            1.0,
            self.y_axis.npts + 1,
        )

        # centered at bins
        self.u_grid = 0.5 * (u_edges[:-1] + u_edges[1:])
        self.v_grid = 0.5 * (v_edges[:-1] + v_edges[1:])

        self.U, self.V = np.meshgrid(self.u_grid, self.v_grid, indexing="xy")
        self.SNG = np.empty_like(self.U, dtype=float)

        # Evaluate the surface through eval(), as requested
        for iy in range(self.y_axis.npts):
            for ix in range(self.x_axis.npts):
                self.SNG[iy, ix] = self.eval_uv_sng_f([self.U[iy, ix], self.V[iy, ix]])

        spmax = np.nanmax(self.SNG)
        self.SNG[self.unsampled_mask] = spmax

    # --------------------------------------------------------------------------

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

        # ----------------------------------------------------------------------
        # Work either in scaled coordinates or in original coordinates
        # ----------------------------------------------------------------------

        data_points = np.array(
            [self.to_scaled([x, y]) for x, y in zip(self.x_data, self.y_data)]
        )

        grid_points = np.column_stack([
            self.x_axis.scale(self.X.ravel()),
            self.y_axis.scale(self.Y.ravel()),
        ])

        dx = 1.0 / self.x_axis.npts
        dy = 1.0 / self.y_axis.npts

        # ----------------------------------------------------------------------
        # Default threshold:
        # a grid point is sampled if there is a data point roughly within
        # one grid-cell diagonal.
        # ----------------------------------------------------------------------

        if max_distance is None:
            max_distance = 0.75 * np.sqrt(dx**2 + dy**2)

        tree = cKDTree(data_points)
        nearest_dist, nearest_idx = tree.query(grid_points, k=1)

        nearest_dist = nearest_dist.reshape(self.X.shape)
        unsampled_mask = nearest_dist > max_distance

        return unsampled_mask, nearest_dist

    # --------------------------------------------------------------------------

    def _build_design_matrix(self, indices):
        """
        Build the RBF design matrix for selected data points.

        Parameters
        ----------
        indices : array-like
            Indices of data points used to build the matrix.

        Returns
        -------
        B : ndarray, shape (ndata_selected, nbasis)
            RBF design matrix.
        """
        indices = np.asarray(indices, dtype=int)

        nxrbf = len(self.x_axis.centers)
        nyrbf = len(self.y_axis.centers)
        nbasis = nxrbf * nyrbf

        B = np.empty((len(indices), nbasis), dtype=float)

        for row, k in enumerate(indices):
            u = self.to_scaled([self.x_data[k], self.y_data[k]])
            B[row, :] = self._basis_2d_uv(u, derivatives=False)

        return B

    # --------------------------------------------------------------------------

    def _solve_svd(self, B, y, rcond):
        """
        Solve a linear least-squares problem using an SVD pseudoinverse.

        Parameters
        ----------
        B : ndarray
            Design matrix.

        y : ndarray
            Target values.

        rcond : float
            Relative singular-value cutoff.

        Returns
        -------
        coeff : ndarray
            Least-squares coefficients.
        """
        U, s, Vt = np.linalg.svd(B, full_matrices=False)

        if len(s) == 0:
            raise RuntimeError("SVD failed: no singular values found")

        cutoff = rcond * np.max(s)

        sinv = np.zeros_like(s)
        mask = s > cutoff
        sinv[mask] = 1.0 / s[mask]

        coeff = Vt.T @ (sinv * (U.T @ y))

        return coeff

    # --------------------------------------------------------------------------

    def fit(self, sx=1.0, sy=1.0, rcond=1.0e-12):
        """
        Fit RBF amplitudes to loaded data using an SVD pseudoinverse.

        This works for both overdetermined and underdetermined systems.
        """
        if self.x_data is None or self.y_data is None or self.e_data is None:
            raise RuntimeError("No data loaded. Call load() first.")

        self.x_axis.width_scale = sx
        self.y_axis.width_scale = sy

        ndata = len(self.e_data)
        nxrbf = len(self.x_axis.centers)
        nyrbf = len(self.y_axis.centers)
        nbasis = nxrbf * nyrbf

        B = np.empty((ndata, nbasis), dtype=float)

        for k in range(ndata):
            u = self.to_scaled([self.x_data[k], self.y_data[k]])
            B[k, :] = self._basis_2d_uv(u, derivatives=False)

        U, s, Vt = np.linalg.svd(B, full_matrices=False)

        smax = np.max(s)
        cutoff = rcond * smax

        sinv = np.zeros_like(s)
        mask = s > cutoff
        sinv[mask] = 1.0 / s[mask]

        coeff = Vt.T @ (sinv * (U.T @ self.e_data))

        self.amplitudes = coeff.reshape((nxrbf, nyrbf))

        # ---------------------------------------------------------------------
        # Calculate fitted values and RMSE on the original data points.
        # ---------------------------------------------------------------------

        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

        return self.amplitudes

    # --------------------------------------------------------------------------

    def fit_width_grid(
        self,
        rcond=1.0e-12,
        width_scales_x=None,
        width_scales_y=None,
        validation_fraction=0.2,
    ):
        """
        Fit RBF amplitudes and optimise separate Gaussian widths for both CVs.

        This scans all combinations of width_scales_x and width_scales_y.
        """
        if self.x_data is None or self.y_data is None or self.e_data is None:
            raise RuntimeError("No data loaded. Call load() first.")

        if width_scales_x is None:
            width_scales_x = np.linspace(0.5, 2.0, 16)

        if width_scales_y is None:
            width_scales_y = np.linspace(0.5, 2.0, 16)

        width_scales_x = np.asarray(width_scales_x, dtype=float)
        width_scales_y = np.asarray(width_scales_y, dtype=float)

        if np.any(width_scales_x <= 0.0) or np.any(width_scales_y <= 0.0):
            raise ValueError("all width scales must be positive")

        ndata = len(self.e_data)


        indices = np.arange(ndata)
        self.rng.shuffle(indices)

        nvalid = max(1, int(round(validation_fraction * ndata)))
        valid_idx = indices[:nvalid]
        train_idx = indices[nvalid:]

        if len(train_idx) == 0:
            raise ValueError("not enough data points for train/validation split")

        best = {
            "sx": None,
            "sy": None,
            "rmse": np.inf,
        }

        for sx in width_scales_x:
            for sy in width_scales_y:
                self.x_axis.width_scale = sx
                self.y_axis.width_scale = sy

                B_train = self._build_design_matrix(train_idx)
                coeff = self._solve_svd(B_train, self.e_data[train_idx], rcond)

                B_valid = self._build_design_matrix(valid_idx)
                e_pred = B_valid @ coeff

                rmse = np.sqrt(np.mean((e_pred - self.e_data[valid_idx]) ** 2))

                print(rmse)

                if rmse < best["rmse"]:
                    best["sx"] = sx
                    best["sy"] = sy
                    best["rmse"] = rmse

        self.x_axis.width_scale = best["sx"]
        self.y_axis.width_scale = best["sy"]

        B = self._build_design_matrix(np.arange(ndata))
        coeff = self._solve_svd(B, self.e_data, rcond)

        nxrbf = len(self.x_axis.centers)
        nyrbf = len(self.y_axis.centers)

        self.amplitudes = coeff.reshape((nxrbf, nyrbf))

        self.best_width_scale_x = best["sx"]
        self.best_width_scale_y = best["sy"]
        self.best_validation_rmse = best["rmse"]

        # ---------------------------------------------------------------------
        # Calculate fitted values and RMSE on the original data points.
        # ---------------------------------------------------------------------

        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

        return self.amplitudes

    # --------------------------------------------------------------------------

    def fit_width_optimize(
        self,
        rcond=1.0e-12,
        width_scales_x=(0.5, 2.5),
        width_scales_y=(0.5, 2.5),
        validation_fraction=0.2,
        verbose=False,
    ):
        """
        Fit RBF amplitudes and optimise separate Gaussian widths for both CVs.

        The width scales are optimised by scipy.optimize.minimize using a
        derivative-free optimiser with bounds.

        Parameters
        ----------
        rcond : float
            Relative cutoff for singular values in the SVD pseudoinverse.

        width_scales_x : tuple(float, float)
            Lower and upper bound for the width scale along CV1.

        width_scales_y : tuple(float, float)
            Lower and upper bound for the width scale along CV2.

        validation_fraction : float
            Fraction of data used for validation.

        random_seed : int
            Seed for reproducible train/validation splitting.

        xtol : float
            Tolerance for width-scale convergence.

        ftol : float
            Tolerance for validation-error convergence.

        maxiter : int
            Maximum number of Powell iterations.

        verbose : bool
            If True, print optimisation progress.

        Returns
        -------
        amplitudes : ndarray
            Fitted RBF amplitudes using the optimised widths and all data.
        """

        if self.x_data is None or self.y_data is None or self.e_data is None:
            raise RuntimeError("No data loaded. Call load() first.")

        if rcond <= 0.0:
            raise ValueError("rcond must be positive")

        if not (0.0 < validation_fraction < 1.0):
            raise ValueError("validation_fraction must be between 0 and 1")

        sx_min, sx_max = map(float, width_scales_x)
        sy_min, sy_max = map(float, width_scales_y)

        if sx_min <= 0.0 or sy_min <= 0.0:
            raise ValueError("width-scale bounds must be positive")

        if sx_max <= sx_min:
            raise ValueError("width_scales_x must be given as (lower, upper)")

        if sy_max <= sy_min:
            raise ValueError("width_scales_y must be given as (lower, upper)")

        ndata = len(self.e_data)

        indices = np.arange(ndata)
        self.rng.shuffle(indices)

        nvalid = max(1, int(round(validation_fraction * ndata)))
        valid_idx = indices[:nvalid]
        train_idx = indices[nvalid:]

        if len(train_idx) == 0:
            raise ValueError("not enough data points for train/validation split")

        def objective(width_scale_vector):
            sx, sy = width_scale_vector

            # Safety guard. Powell should respect bounds, but this keeps
            # the objective well-defined even if trial points drift slightly.
            if sx <= 0.0 or sy <= 0.0:
                return np.inf

            self.x_axis.width_scale = sx
            self.y_axis.width_scale = sy

            try:
                B_train = self._build_design_matrix(train_idx)
                coeff = self._solve_svd(B_train, self.e_data[train_idx], rcond)

                B_valid = self._build_design_matrix(valid_idx)
                e_pred = B_valid @ coeff

                rmse = np.sqrt(np.mean((e_pred - self.e_data[valid_idx]) ** 2))

            except Exception:
                rmse = np.inf

            print(f"  sx = {sx:12.6f}  sy = {sy:12.6f}  RMSE = {rmse:12.6f}")

            return rmse

        x0 = np.array(
            [
                0.5 * (sx_min + sx_max),
                0.5 * (sy_min + sy_max),
            ],
            dtype=float,
        )

        bounds = [
            (sx_min, sx_max),
            (sy_min, sy_max),
        ]

        result = minimize(
            objective,
            x0,
            method="Nelder-Mead",
            bounds=bounds,
        )

        if not result.success and verbose:
            print("Width optimisation warning:", result.message)

        best_sx, best_sy = result.x

        self.x_axis.width_scale = best_sx
        self.y_axis.width_scale = best_sy

        # Final refit using all data and the optimised widths.
        B = self._build_design_matrix(np.arange(ndata))
        coeff = self._solve_svd(B, self.e_data, rcond)

        nxrbf = len(self.x_axis.centers)
        nyrbf = len(self.y_axis.centers)

        self.amplitudes = coeff.reshape((nxrbf, nyrbf))

        self.best_width_scale_x = float(best_sx)
        self.best_width_scale_y = float(best_sy)
        self.best_validation_rmse = float(result.fun)
        self.width_optimization_result = result

        # ---------------------------------------------------------------------
        # Calculate fitted values and RMSE on the original data points.
        # ---------------------------------------------------------------------

        self.e_fit = B @ coeff
        self.fit_residuals = self.e_fit - self.e_data
        self.fit_rmse = float(np.sqrt(np.mean(self.fit_residuals**2)))

        return self.amplitudes

    # --------------------------------------------------------------------------
    # RBF basis
    # --------------------------------------------------------------------------

    def _basis_1d(self, u, axis):
        """
        Return 1D Gaussian basis values and derivatives in scaled coordinates.

        Returns
        -------
        g : ndarray
            Gaussian values.
        dg : ndarray
            First derivatives with respect to scaled coordinate.
        d2g : ndarray
            Second derivatives with respect to scaled coordinate.
        """
        du = axis.delta(u)
        s = axis.width * axis.width_scale

        g = np.exp(-0.5 * (du / s) ** 2)

        dg = g * (-du / s**2)

        d2g = g * ((du**2 / s**4) - (1.0 / s**2))

        return g, dg, d2g

    # --------------------------------------------------------------------------

    def _basis_2d_uv(self, u, derivatives=False):
        """
        Evaluate all 2D RBFs in scaled coordinates.

        Parameters
        ----------
        u : array-like, shape (2,)
            Scaled coordinates [u, v].

        derivatives : bool
            If True, also return first and second derivatives in scaled coordinates.

        Returns
        -------
        If derivatives is False:

            phi

        If derivatives is True:

            phi, phi_u, phi_v, phi_uu, phi_vv, phi_uv
        """
        u = np.asarray(u, dtype=float)

        if u.shape != (2,):
            raise ValueError("u must be a vector with two elements: [u, v]")

        ux = u[0]
        uy = u[1]

        gx, dgx, d2gx = self._basis_1d(ux, self.x_axis)
        gy, dgy, d2gy = self._basis_1d(uy, self.y_axis)

        phi = np.outer(gx, gy)

        if not derivatives:
            return phi.ravel()

        phi_u = np.outer(dgx, gy)
        phi_v = np.outer(gx, dgy)

        phi_uu = np.outer(d2gx, gy)
        phi_vv = np.outer(gx, d2gy)
        phi_uv = np.outer(dgx, dgy)

        return (
            phi.ravel(),
            phi_u.ravel(),
            phi_v.ravel(),
            phi_uu.ravel(),
            phi_vv.ravel(),
            phi_uv.ravel(),
        )

# ==============================================================================
# Evaluation
# ==============================================================================

    def eval_uv(self, u):
        """
        Evaluate the surface in scaled coordinates.

        Parameters
        ----------
        u : array-like, shape (2,)
            Scaled coordinates [u, v].

        Returns
        -------
        value : float
            Energy E(u, v).

        gradient : ndarray, shape (2,)
            Gradient in scaled coordinates:

                [dE/du, dE/dv]

        hessian : ndarray, shape (2, 2)
            Hessian in scaled coordinates:

                [[d2E/du2,  d2E/dudv],
                [d2E/dvdu, d2E/dv2 ]]
        """
        if self.amplitudes is None:
            raise RuntimeError("The surface has not been fitted. Call fit() first.")

        A = self.amplitudes.ravel()

        phi, phi_u, phi_v, phi_uu, phi_vv, phi_uv = self._basis_2d_uv(
            u, derivatives=True
        )

        value = float(np.dot(A, phi))

        gradient = np.array(
            [
                np.dot(A, phi_u),
                np.dot(A, phi_v),
            ],
            dtype=float,
        )

        hessian = np.array(
            [
                [np.dot(A, phi_uu), np.dot(A, phi_uv)],
                [np.dot(A, phi_uv), np.dot(A, phi_vv)],
            ],
            dtype=float,
        )

        return value, gradient, hessian

    # --------------------------------------------------------------------------

    def eval(self, x):
        """
        Evaluate the surface in original coordinates.

        Parameters
        ----------
        x : array-like, shape (2,)
            Original coordinates [x, y].

        Returns
        -------
        value : float
            Energy E(x, y).

        gradient : ndarray, shape (2,)
            Gradient in original coordinates:

                [dE/dx, dE/dy]

        hessian : ndarray, shape (2, 2)
            Hessian in original coordinates:

                [[d2E/dx2,  d2E/dxdy],
                [d2E/dydx, d2E/dy2 ]]
        """
        x = np.asarray(x, dtype=float)

        if x.shape != (2,):
            raise ValueError("x must be a vector with two elements: [x, y]")

        u = self.to_scaled(x)

        value, gradient_uv, hessian_uv = self.eval_uv(u)

        rx = self.x_axis.range
        ry = self.y_axis.range

        scale = np.array([rx, ry], dtype=float)

        gradient = gradient_uv / scale

        hessian = hessian_uv / np.outer(scale, scale)

        return value, gradient, hessian

    # --------------------------------------------------------------------------

    def eval_uv_sng_f(self, u):

        """
        Evaluate the squared norm of energy gradient (SNG). SNG is zero at all stationary points and positive elsewhere.

        Returns
        -------
            value of SNG 
        """

        value, gradient, hessian = self.eval_uv(u)

        sneg_v = gradient[0]**2 + gradient[1]**2

        return sneg_v
    
    # --------------------------------------------------------------------------

    def eval_uv_sng_fg(self, u):

        """
        Evaluate the squared norm of energy gradient (SNG). SNG is zero at all stationary points and positive elsewhere.

        Returns
        -------
         value of SNG and its gradient
        """

        value, gradient, hessian = self.eval_uv(u)

        sneg_v = gradient[0]**2 + gradient[1]**2

        sneg_g = np.array(
            [
                2.0 * gradient[0] * hessian[0][0] + 2.0 * gradient[0] * hessian[0][1],
                2.0 * gradient[1] * hessian[1][0] + 2.0 * gradient[1] * hessian[1][1]
            ],
            dtype=float,
        )

        return sneg_v, sneg_g

# ==============================================================================
# Plots
# ==============================================================================

    def make_gnuplot_like_colormaps(self,zmin,zmax,contour_spacing):
        """
        Create matplotlib colormaps compatible with the gnuplot palette:

            set palette maxcolors 18
            set palette defined (
                0 '#3288BD',
                1 '#66C2A5',
                2 '#ABDDA4',
                3 '#E6F598',
                4 '#FEE08B',
                5 '#FDAE61',
                6 '#F46D43',
                7 '#D53E4F'
            )

        Returns
        -------
        fes_cmap : Colormap
            Discrete color map for FES.
        fes_norm : BoundaryNorm
            contour_spacing binning from zmin to zmax.
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

        # contour_spacing bins
        bounds = np.arange(zmin, zmax + contour_spacing, contour_spacing)

        N = len(bounds)

        fes_norm = BoundaryNorm(bounds, N, clip=True)

        # Interpolate the 8 anchor colors to 18 discrete colors
        fes_cmap = LinearSegmentedColormap.from_list("gnuplot_like_fes", gnuplot_colors, N)

        # Watershed regions are categorical, so use the anchor colors directly
        ws_cmap = ListedColormap(gnuplot_colors, name="gnuplot_like_regions")

        return fes_cmap, fes_norm, ws_cmap
    
    # --------------------------------------------------------------------------

    def _finish_plot(self, fig, ax, show, save, dpi):

        ax.set_xlabel(self.x_axis.label)
        ax.set_ylabel(self.y_axis.label)

        ax.set_xlim(self.x_axis.cvmin, self.x_axis.cvmax)
        ax.set_ylim(self.y_axis.cvmin, self.y_axis.cvmax)

        ax.set_aspect('equal', adjustable='box')
        fig.tight_layout()

        if save is not None:
            fig.savefig(save, dpi=dpi)
        if show:
            plt.show()

        plt.close()

    # --------------------------------------------------------------------------

    def plot_raw(
        self,
        marker_size=25,
        figsize=None,
        show=True,
        save=None,
        dpi=300,
    ):
        """
        Visualise the original 2D energy surface.

        Parameters
        ----------
        show : bool
            If True, call plt.show().

        save : str or None
            If not None, save the figure to this filename.
        """

        if self.e_data is None:
            raise RuntimeError("The surface has not been loaded. Call load() first.")

        fes_cmap, fes_norm, ws_cmap = self.make_gnuplot_like_colormaps(0.0,self.zmax,self.contour_spacing)

        fig, ax = plt.subplots(figsize=figsize)

        x = np.asarray(self.x_data, dtype=float)
        y = np.asarray(self.y_data, dtype=float)
        e = np.asarray(self.e_data, dtype=float)

        if self.zmax is not None:
            eplot = np.minimum(e, self.zmax)
        else:
            eplot = e

        # plot the FES
        sc = ax.scatter(
            x,
            y,
            c=eplot,
            s=marker_size,
            cmap=fes_cmap,
            norm=fes_norm,
            edgecolors="none",
        )

        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(self.ene_label)

        ax.set_title("Original FES")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_rbf(
        self,
        figsize=None,
        show=True,
        save=None,
        dpi=300,
    ):
        """
        Visualise the fitted 2D energy surface using eval().

        Parameters
        ----------
        show : bool
            If True, call plt.show().

        save : str or None
            If not None, save the figure to this filename.
        """

        if self.amplitudes is None:
            raise RuntimeError("The surface has not been fitted. Call fit() first.")

        fes_cmap, fes_norm, ws_cmap = self.make_gnuplot_like_colormaps(0.0,self.zmax,self.contour_spacing)

        fig, ax = plt.subplots(figsize=figsize)

        # plot the FES
        im = ax.pcolormesh(
            self.X,
            self.Y,
            self.Z,
            shading="auto",
            cmap=fes_cmap,
            norm=fes_norm,
        )

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(self.ene_label)

        # Contours
        cmin = 0.0
        cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing

        levels = np.arange(cmin, cmax + self.contour_spacing, self.contour_spacing)

        contours = ax.contour(
            self.X,
            self.Y,
            self.Z,
            levels=levels,
            colors="black",
            linewidths=0.6,
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt="%.0f",
        )

        ax.set_title("Interpolated FES")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_uv_sng(
        self,
        zmax=None,
        contour_spacing=1.0,
        cmap="viridis",
        figsize=None,
        show=True,
        save=None,
        dpi=300,
    ):
        """
        Visualise the fitted 2D energy surface using eval().

        Parameters
        ----------
        contour_spacing : float
            Spacing between contour lines in energy units.

        zmax : float or None
            Maximum energy shown in the map and contours.
            If None, the full energy range is used.

        cmap : str
            Matplotlib colormap name.

        show : bool
            If True, call plt.show().

        save : str or None
            If not None, save the figure to this filename.
        """

        if self.SNG is None:
            raise RuntimeError("The SNG surface has not been calculated. Call calc_ene_and_sng() first.")

        if contour_spacing <= 0.0:
            raise ValueError("contour_spacing must be positive")

        if zmax is not None:
            zmax = float(zmax)
            Zplot = np.minimum(self.SNG, zmax)
        else:
            Zplot = self.SNG

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.pcolormesh(
            self.U,
            self.V,
            Zplot,
            shading="auto",
            cmap=cmap,
            vmin=np.nanmin(Zplot),
            vmax=np.nanmax(Zplot),
        )

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("SNG Metric")

        zmin_plot = np.nanmin(Zplot)
        zmax_plot = np.nanmax(Zplot)

        cmin = np.floor(zmin_plot / contour_spacing) * contour_spacing
        cmax = np.ceil(zmax_plot / contour_spacing) * contour_spacing

        levels = np.arange(cmin, cmax + contour_spacing, contour_spacing)

        contours = ax.contour(
            self.X,
            self.Y,
            Zplot,
            levels=levels,
            colors="black",
            linewidths=0.6,
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt="%.1f",
        )

        ax.set_xlabel("CV1 [scaled]")
        ax.set_ylabel("CV2 [scaled]")
        ax.set_title("SNG Surface")

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)

        fig.tight_layout()

        if save is not None:
            fig.savefig(save, dpi=dpi)

        if show:
            plt.show()

        plt.close()

    # --------------------------------------------------------------------------

    def plot_rbf_with_sp_guesses(
        self,
        figsize=None,
        show=True,
        save=None,
        dpi=300,
    ):
        """
        Visualise the fitted 2D energy surface using eval() + sp guesses.

        Parameters
        ----------
        show : bool
            If True, call plt.show().

        save : str or None
            If not None, save the figure to this filename.
        """

        if self.amplitudes is None:
            raise RuntimeError("The surface has not been fitted. Call fit() first.")

        fes_cmap, fes_norm, ws_cmap = self.make_gnuplot_like_colormaps(0.0,self.zmax,self.contour_spacing)

        fig, ax = plt.subplots(figsize=figsize)

        # plot the FES
        im = ax.pcolormesh(
            self.X,
            self.Y,
            self.Z,
            shading="auto",
            cmap=fes_cmap,
            norm=fes_norm,
        )

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(self.ene_label)

        # Contours
        cmin = 0.0
        cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing

        levels = np.arange(cmin, cmax + self.contour_spacing, self.contour_spacing)

        contours = ax.contour(
            self.X,
            self.Y,
            self.Z,
            levels=levels,
            colors="black",
            linewidths=0.6,
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt="%.0f",
        )

        # -------------------------------------------------------------------------
        # Overlay detected zero/stationary regions
        # -------------------------------------------------------------------------

        if len(self.sp_guesses) > 0:

            rx = []
            ry = []
            labels = []

            for k, region in enumerate(self.sp_guesses):
                rx.append(region['x'])
                ry.append(region['y'])
                labels.append(region['id'])

            ax.scatter(
                rx,
                ry,
                s=60,
                marker="o",
                facecolors="none",
                edgecolors="red",
                linewidths=1.8,
                label="detected stationary points",
                zorder=10,
            )

            ax.scatter(
                rx,
                ry,
                s=15,
                marker='x',
                color="red",
                zorder=11,
            )

            for k, (px, py) in enumerate(zip(rx, ry)):
                label = str(labels[k])

                ax.text(
                    px,
                    py,
                    " " + label,
                    color="red",
                    fontsize=9,
                    ha="left",
                    va="center",
                    zorder=12,
                )

            ax.legend(loc="best")

        # -------------------------------------------------------------------------

        ax.set_title("Interpolated FES with Stationary Point Guesses")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_rbf_with_sp_optimized(
        self,
        figsize=None,
        show=True,
        save=None,
        dpi=300,
        point_size=70,
        ellipse_lw=2.0,
        ellipse_npts=181,
        point_color="red",
        ellipse_color="red",
        label_color="red",
        label_offset=(0.05, 0.05),
    ):
        """
        Plot interpolated FES and stationary points with labelled ellipses.

        Parameters
        ----------
        show : bool
            If True, call plt.show().

        save : str or None
            If not None, save the figure to this filename.

        point_size : float
            Marker size for stationary points.

        ellipse_lw : float
            Line width for ellipses.

        ellipse_npts : int
            Number of points used to draw each ellipse.

        point_color, ellipse_color, label_color : str
            Matplotlib colors.

        label_offset : tuple(float, float)
            Offset added to labels in original-coordinate units.
        """

        if self.amplitudes is None:
            raise RuntimeError("The surface has not been fitted. Call fit() first.")

        if self.X is None or self.Y is None or self.Z is None:
            raise RuntimeError("Interpolated FES is not available. Call calc_ene_and_sng() first.")

        if ellipse_npts < 4:
            raise ValueError("ellipse_npts must be at least 4.")

        # ---------------------------------------------------------------------
        # Prepare plot.
        # ---------------------------------------------------------------------

        fes_cmap, fes_norm, ws_cmap = self.make_gnuplot_like_colormaps(
            0.0,
            self.zmax,
            self.contour_spacing,
        )

        fig, ax = plt.subplots(figsize=figsize)

        # ---------------------------------------------------------------------
        # Plot interpolated FES.
        # ---------------------------------------------------------------------

        im = ax.pcolormesh(
            self.X,
            self.Y,
            self.Z,
            shading="auto",
            cmap=fes_cmap,
            norm=fes_norm,
        )

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(self.ene_label)

        cmin = 0.0
        cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing
        levels = np.arange(cmin, cmax + self.contour_spacing, self.contour_spacing)

        contours = ax.contour(
            self.X,
            self.Y,
            self.Z,
            levels=levels,
            colors="black",
            linewidths=0.6,
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt="%.0f",
        )

        # ---------------------------------------------------------------------
        # Overlay stationary points and ellipses.
        # ---------------------------------------------------------------------

        t_values = np.linspace(0.0, 2.0 * np.pi, ellipse_npts)

        for p in self.sp_optimized:
            x0 = p['x']
            y0 = p['y']
            r1 = p['r1']
            r2 = p['r2']
            angle = p['angle'] / rad2deg

            ca = np.cos(angle)
            sa = np.sin(angle)

            xe = x0 + r1 * np.cos(t_values) * ca - r2 * np.sin(t_values) * sa
            ye = y0 + r1 * np.cos(t_values) * sa + r2 * np.sin(t_values) * ca

            ax.plot(
                xe,
                ye,
                color=ellipse_color,
                linewidth=ellipse_lw,
                zorder=10,
            )

            ax.scatter(
                [x0],
                [y0],
                s=point_size,
                marker="o",
                facecolors="none",
                edgecolors=point_color,
                linewidths=1.8,
                zorder=11,
            )

            ax.scatter(
                [x0],
                [y0],
                s=20,
                marker='x',
                color=point_color,
                zorder=12,
            )

            label = str(p['id'])

            if p['type'] != "":
                label = f"{label}:{p['type']}"

            ax.text(
                x0 + label_offset[0],
                y0 + label_offset[1],
                label,
                color=label_color,
                fontsize=9,
                ha="left",
                va="bottom",
                zorder=13,
            )

        # ---------------------------------------------------------------------
        # Final styling.
        # ---------------------------------------------------------------------

        ax.set_title("Interpolated FES with Optimized Stationary Points")
        self._finish_plot(fig, ax, show, save, dpi)

    # --------------------------------------------------------------------------

    def plot_rbf_with_watershed_basins(
        self,
        figsize=None,
        show=True,
        save=None,
        dpi=300,
    ):

        # ---------------------------------------------------------------------
        # Prepare plot.
        # ---------------------------------------------------------------------

        fes_cmap, fes_norm, ws_cmap = self.make_gnuplot_like_colormaps(
            0.0,
            self.zmax,
            self.contour_spacing,
        )

        fig, ax = plt.subplots(figsize=figsize)

        # ---------------------------------------------------------------------
        # Plot interpolated FES.
        # ---------------------------------------------------------------------

        im = ax.pcolormesh(
            self.X,
            self.Y,
            self.Z,
            shading="auto",
            cmap=fes_cmap,
            norm=fes_norm,
            zorder=1,
        )

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(self.ene_label)

        cmin = 0.0
        cmax = np.ceil(self.zmax / self.contour_spacing) * self.contour_spacing
        levels = np.arange(cmin, cmax + self.contour_spacing, self.contour_spacing)

        contours = ax.contour(
            self.X,
            self.Y,
            self.Z,
            levels=levels,
            colors="black",
            linewidths=0.6,
            zorder=2,
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=8,
            fmt="%.0f",
            zorder=3,
        )

        # -------------------------------------------------------------------------
        # Watershed regions (categorical overlay)
        # -------------------------------------------------------------------------

        # Map region indices cyclically to the 8 palette anchor colors
        region_img = np.zeros_like(self.basins_labels, dtype=float)
        mask = self.basins_labels > 0
        region_img[mask] = ((self.basins_labels[mask] - 1) % ws_cmap.N) + 1
        region_img = np.ma.masked_where(self.basins_labels == 0, region_img)

        ax.imshow(
            region_img,
            origin="lower",
            extent=(self.x_axis.cvmin, self.x_axis.cvmax, self.y_axis.cvmin, self.y_axis.cvmax),
            interpolation="nearest",
            cmap=ws_cmap,
            alpha=0.80,
            aspect="auto",
            vmin=1,
            vmax=ws_cmap.N,
            zorder=4,
        )

        # -------------------------------------------------------------------------
        # Watershed boundaries
        # -------------------------------------------------------------------------

        if np.max(self.basins_labels) > 0:
            ax.contour(
                self.X,
                self.Y,
                self.basins_labels,
                levels=np.arange(0.5, self.basins_labels.max() + 0.5, 1.0),
                colors="white",
                linewidths=1.0,
                alpha=0.90,
                zorder=5,
            )

        # -------------------------------------------------------------------------
        # Local-minimum seeds
        # -------------------------------------------------------------------------

        for i, p in enumerate(self.sp_optimized, start=1):
            ax.plot(p['x'], p['y'], "ko", markersize=4,zorder=6)
            ax.text(
                p['x'],
                p['y'],
                f" {i}",
                color="white",
                fontsize=8,
                weight="bold",
                ha="left",
                va="center",
                zorder=6,
            )

        # ---------------------------------------------------------------------
        # Final styling.
        # ---------------------------------------------------------------------

        ax.set_title("Interpolated FES with Minima Basins")
        self._finish_plot(fig, ax, show, save, dpi)

# ==============================================================================
# Detect stationary points
# ==============================================================================

    def detect_sp_regions_connected(
        self,
        eps=10,
        connectivity=2,
        representative="minabs",
        min_size=1,
    ):
        """
        Detect connected regions of nearly-zero values on a regular 2D grid.

        Parameters
        ----------
        eps : float
            Threshold for nearly-zero values.
        connectivity : int
            1 = 4-neighbour connectivity
            2 = 8-neighbour connectivity
        representative : str
            "minabs"  : grid point with smallest |F|
            "centroid": grid point closest to region centroid
        min_size : int
            Ignore regions with fewer than min_size grid points.

        Returns
        -------
        self.sp_guesses : list of dict
            One dictionary per detected SP.
        self.sp_guess_labels : 2D array
            Integer region labels. Background is 0.
        sp_guess_mask : 2D bool array
            Thresholded nearly-zero mask.
        """

        self.sp_guess_mask = (self.SNG >= 0.0) & (self.SNG <= eps) & (self.Z < self.zmax)

        if connectivity == 1:
            structure = ndimage.generate_binary_structure(2, 1)  # 4-neighbour
        elif connectivity == 2:
            structure = ndimage.generate_binary_structure(2, 2)  # 8-neighbour
        else:
            raise ValueError("connectivity must be 1 or 2")

        self.sp_guess_labels, nlabels = ndimage.label(self.sp_guess_mask, structure=structure)

        self.sp_guesses = []

        for label in range(1, nlabels + 1):
            idx = np.argwhere(self.sp_guess_labels == label)

            if idx.shape[0] < min_size:
                self.sp_guess_labels[self.sp_guess_labels == label] = 0
                continue

            values = self.SNG[self.sp_guess_labels == label]

            if representative == "minabs":
                k = np.argmin(np.abs(values))
                i, j = idx[k]

            elif representative == "centroid":
                centroid = idx.mean(axis=0)
                distances2 = np.sum((idx - centroid) ** 2, axis=1)
                k = np.argmin(distances2)
                i, j = idx[k]

            else:
                raise ValueError("representative must be 'minabs' or 'centroid'")

            self.sp_guesses.append({
                'id': label,
                "i": int(i),
                "j": int(j),
                "u": float(self.U[i, j]),
                "v": float(self.V[i, j]),
                'sng': float(self.SNG[i, j]),
                'x': float(self.X[i, j]),
                'y': float(self.Y[i, j]),
                'ene': float(self.Z[i, j]),
                "size": int(idx.shape[0]),
            })

        return self.sp_guesses, self.sp_guess_labels, self.sp_guess_mask

# ==============================================================================
# Optimize stationary points
# ==============================================================================

    def find_sp(self, gx, gy, trust_r0=0.06):
        """
        Locate SP on FES.

        Parameters
        ----------
        gx, gy : str or float
            Initial guess location of stationary point.
        """

        # -------------------------------------------------------------------------
        # Input guess
        # -------------------------------------------------------------------------

        x0s = np.array([float(gx), float(gy)], dtype=float)
        x0s = self.to_scaled(x0s)

        # -------------------------------------------------------------------------
        # Find SP
        # -------------------------------------------------------------------------

        best_res = None
        best_fun = np.inf
        best_x0 = None

        print("")
        print("Stationary point scan:")

        # we need to search in a small box in the case of very narrow minima
        bounds = [(x0s[0]-trust_r0, x0s[0]+trust_r0), (x0s[1]-trust_r0, x0s[1]+trust_r0)]

        # try gradient based optimizer
        for attempt in range(0, 10):

            if attempt > 0:
                dx = self.rng.uniform(
                    low=-trust_r0/3.0,  # use smaller interval
                    high=trust_r0/3.0,
                    size=x0s.shape
                )

                x0 = np.zeros(2)
                x0 = x0s + dx
            else:
                x0 = x0s

            print(f"  x0 = {x0}")

            res = minimize(
                self.eval_uv_sng_fg,
                x0,
                method="L-BFGS-B",
                jac=True,
            )

            print(f"      err = {res.fun:12.6e}, x = {res.x}")

            if (res.x[0] < 0.0) or (res.x[0] > 1.0) or (res.x[1] < 0.0) or (res.x[1] > 1.0):
                print(" >> solution out-of-box: ignoring")
                continue

            if np.linalg.norm(res.x-x0) > trust_r0:
                print(" >> solution out-of-trust-region: ignoring")
                continue 

            if res.fun < best_fun:
                best_fun = res.fun
                best_res = res
                best_x0 = x0.copy()

        if best_res is None:
            # try non-gradient optimizers
            for attempt in range(0, 10):

                if attempt > 0:
                    dx = self.rng.uniform(
                        low=-trust_r0/3.0,  # use smaller interval
                        high=trust_r0/3.0,
                        size=x0s.shape
                    )

                    x0 = np.zeros(2)
                    x0 = x0s + dx
                else:
                    x0 = x0s

                print(f"  x0 = {x0}")

                res = minimize(
                    self.eval_uv_sng_f,
                    x0,
                    method="Nelder-Mead",
                    bounds=bounds
                )

                print(f"      err = {res.fun:12.6e}, x = {res.x}")

                if (res.x[0] < 0.0) or (res.x[0] > 1.0) or (res.x[1] < 0.0) or (res.x[1] > 1.0):
                    print(" >> solution out-of-box: ignoring")
                    continue

                if np.linalg.norm(res.x-x0) > trust_r0:
                    print(" >> solution out-of-trust-region: ignoring")
                    continue 

                if res.fun < best_fun:
                    best_fun = res.fun
                    best_res = res
                    best_x0 = x0.copy()


        if best_res is None:
            # not found
            item = {
                'id': -1
            }
            return item

        print("")
        print("Optimized SP parameters:")
        print(f"  initial x0      = {best_x0}")
        print(f"  best error      = {best_fun:12.6e}")
        print(f"  best parameters = {best_res.x}")

        print("")
        print(best_res)

        self.xopt = best_res.x

        print("")
        print(f"xopt = {self.xopt}")
        print(f"scaled xopt = {self.to_scaled(self.xopt)}")

        # -------------------------------------------------------------------------
        # Analyze shape of stationary point from Hessian at xopt
        # -------------------------------------------------------------------------

        [self.xoptene, g, hess] = self.eval_uv(self.xopt)

        eigvals, eigvecs = np.linalg.eigh(hess)

        print("")
        print("Hessian eigenvalues:\n",eigvals)

        idx = np.argsort(eigvals)
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        print("")
        num_negative = np.sum(eigvals < 0.0)
        if num_negative == 2:
            print("Type of SP: local maximum")
            xopttype = "M"
            self.ts_ellipse = 0
        elif num_negative == 1:
            print("Type of SP: transition state")
            xopttype = "T"
            self.ts_ellipse = 1
        else:
            print("Type of SP: local minimum")
            xopttype = "S"
            self.ts_ellipse = 0

        # -------------------------------------------------------------------------
        # Optimize size and rotation of ellipse
        # -------------------------------------------------------------------------

        # # initial params
        # r1 = np.sqrt(2.0 * self.thr / abs(eigvals[0]))
        # r2 = np.sqrt(2.0 * self.thr / abs(eigvals[1]))
        #
        # # Orientation of the first eigenvector
        # angle = np.arctan2(eigvecs[1, 0], eigvecs[0, 0])

        r1 = 0.01
        r2 = 0.02
        angle = 0.0

        best_res = None
        best_fun = np.inf
        best_x0 = None

        print("")
        print("Ellipse scan:")

        for angle_deg in range(0, 125, 30):
            angle_rad = angle + angle_deg * np.pi / 180.0

            x0_ellipse = np.array([r1, r2, angle_rad], dtype=float)

            print(f"  angle = {angle_deg:6.1f} deg, x0 = {x0_ellipse}")

            # first rough optimization
            res_nm = minimize(
                self.my_ellipse_err,
                x0_ellipse,
                method="Nelder-Mead"
            )

            print(f"      err = {res_nm.fun:12.6e}, x = {res_nm.x}")

            if res_nm.fun < best_fun:
                best_fun = res_nm.fun
                best_res = res_nm
                best_x0 = x0_ellipse.copy()

        print("")
        print(best_res)

        print("")
        print("Optimized ellipse parameters:")
        print(f"  initial x0      = {best_x0}")
        print(f"  best error      = {best_fun:12.6e}")
        print(f"  best parameters = {best_res.x}")

        sr1 = best_res.x[0]
        sr2 = best_res.x[1]
        angle = best_res.x[2]

        sng = self.eval_uv_sng_f(self.xopt)

        xopt, r1, r2, angle, eigvecs_orig = self.ellipse_scaled_to_original(self.xopt, sr1, sr2, angle)

        sp_label = len(self.sp_optimized) + 1

        item = {
            'id': sp_label,
            'x': float(xopt[0]),
            'y': float(xopt[1]),
            'ene': float(self.xoptene),
            'sng': float(sng),
            'r1': float(r1),
            'r2': float(r2),
            'sr1': float(sr1),
            'sr2': float(sr2),
            'angle': float(angle * rad2deg),
            'type': xopttype,
            'sl1': eigvals[0],
            'sl2': eigvals[1],
        }

        self.sp_optimized.append(item)

        return item

    # --------------------------------------------------------------------------

    def remove_near_duplicate_sp_optimized(self, min_distance_uv=0.01):
        """
        Remove nearly duplicate optimized stationary points.

        Two points are considered duplicates only if they have the same type
        ('S', 'T', or 'M') and their distance in scaled UV coordinates is
        smaller than or equal to min_distance_uv.

        From each duplicate group, the point with the lowest energy is kept.

        Parameters
        ----------
        min_distance_uv : float
            Minimum distance between two SPs to be considered as individual points.

        Returns
        -------
        removed : list of dict
            List of removed stationary points.
        """

        print("")
        print("# Detecting near duplicate points ...")

        if min_distance_uv <= 0.0:
            raise ValueError("min_distance_uv must be positive")

        if len(self.sp_optimized) <= 1:
            return []

        # ------------------------------------------------------------------
        # Sort points by increasing energy.
        # Therefore, when a conflict is found, the already kept point always
        # has lower or equal energy.
        # ------------------------------------------------------------------

        kept = []
        removed = []

        for pidx, p in enumerate(self.sp_optimized):

            p_type = p.get('type', "")
            if p_type not in ("S", "T", "M"):
                # Unknown type: do not compare it with S/T/M points.
                kept.append(p)
                continue

            p_uv = self.to_scaled([p['x'], p['y']])

            keep_item = p

            if keep_item in kept:
                continue

            duplicate_of = None

            for q in self.sp_optimized[pidx+1:-1]:

                q_type = q.get('type', "")
                if q_type != p_type:
                    continue

                q_uv = self.to_scaled([q['x'], q['y']])

                du = p_uv - q_uv

                dist_uv = np.linalg.norm(du)

                if dist_uv <= min_distance_uv:
                    if p["type"] == "S":
                        if p["ene"] > q["ene"]:
                            keep_item = q
                            duplicate_of = p
                            break
                    if p["type"] == "T":
                        if p["ene"] < q["ene"]:
                            keep_item = q
                            duplicate_of = p
                            break

            if duplicate_of is None:
                kept.append(keep_item)
            else:
                kept.append(keep_item)
                removed.append(duplicate_of)
                print(
                    "  Removing duplicate SP: "
                    f"label = {duplicate_of['id']}, "
                    f"type = {duplicate_of['type']}, "
                    f"ene = {duplicate_of['ene']:.6f}, "
                    f"sng = {duplicate_of['sng']:.6f}, "
                    f"x = {duplicate_of['x']:.6f}, "
                    f"y = {duplicate_of['y']:.6f}; "
                    f"kept label = {keep_item['id']}, "
                    f"kept ene = {keep_item['ene']:.6f}"
                )

        # ------------------------------------------------------------------
        # Restore a stable order, preferably by label.
        # ------------------------------------------------------------------

        kept = sorted(kept, key=lambda p: int(p['id']))

        self.sp_optimized = kept

        print(f"  Removed {len(removed)} near-duplicate optimized stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return removed
    
    # --------------------------------------------------------------------------

    def remove_stationary_point_outliers(self, minslam=5.0, maxslam=5000.0):

        print(f"")
        print(f"# Detecting Hessian outliers ...")
        print(f"  Allowed Hessian eigenvalues: <{minslam:.1f}, {maxslam:.1f}>")
    
        kept = []
        removed = []

        for item in self.sp_optimized:
            outlier = (
                min(abs(item["sl1"]), abs(item["sl2"])) < minslam
                or max(abs(item["sl1"]), abs(item["sl2"])) > maxslam
            )
            if outlier:
                print(f"  Removing stationary point {item['id']} ({item['type']}) due to Hessian eigenvalues ({item['sl1']:.1f},{item['sl2']:.1f}) out-of-allowed values.")
                removed.append(item)
            else:
                kept.append(item)

        self.sp_optimized = kept

        print(f"  Removed {len(removed)} stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return
    
# --------------------------------------------------------------------------

    def remove_small_ellipses(self, minr=0.005):

        print(f"")
        print(f"# Detecting small ellipses ...")
        print(f"  Minimum radius: {minr:.4f}")
    
        kept = []
        removed = []

        for item in self.sp_optimized:
            outlier = (
                min(abs(item["sr1"]), abs(item["sr2"])) < minr
            )
            if outlier:
                print(f"  Removing stationary point {item['id']} ({item['type']}) due to small ellipses ({item['sr1']:.4f},{item['sr2']:.4f}) out-of-allowed values.")
                removed.append(item)
            else:
                kept.append(item)

        self.sp_optimized = kept

        print(f"  Removed {len(removed)} stationary point(s).")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return
    
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

        return

    # --------------------------------------------------------------------------

    def ellipse_scaled_to_original(self, xopt, r1, r2, angle):
        """
        Convert ellipse parameters from scaled coordinates to original coordinates.

        Parameters
        ----------
        xopt : array-like
            Ellipse center in scaled coordinates [u1, u2].
        r1, r2 : float
            Semi-axis lengths in scaled coordinates.
        angle : float
            Rotation angle in scaled coordinates, in radians.

        Returns
        -------
        xopt_orig : ndarray
            Ellipse center in original coordinates.
        R1_orig, R2_orig : float
            Semi-axis lengths in original coordinates.
        angle_orig : float
            Rotation angle of the first principal axis in original coordinates, radians.
        eigvecs_orig : ndarray
            Eigenvectors / principal directions in original coordinates.
        """

        # coordinate scaling factors
        sx = self.x_axis.cvmax - self.x_axis.cvmin
        sy = self.y_axis.cvmax - self.y_axis.cvmin

        # center
        xopt_orig = self.from_scaled(np.asarray(xopt, dtype=float))

        # rotation matrix in scaled coordinates
        ca = np.cos(angle)
        sa = np.sin(angle)

        R = np.array([
            [ca, -sa],
            [sa,  ca]
        ])

        # diagonal matrix of scaled radii
        D = np.diag([r1, r2])

        # coordinate scaling matrix: original = offset + S * scaled
        S = np.diag([sx, sy])

        # Linear map from unit circle to original-coordinate ellipse
        #
        # scaled ellipse:
        #     u = xopt + R @ D @ q
        #
        # original ellipse:
        #     x = xopt_orig + S @ R @ D @ q
        #
        M = S @ R @ D

        # Shape matrix in original coordinates
        C = M @ M.T

        # Eigenvalues of C are squared semi-axis lengths
        eigvals, eigvecs = np.linalg.eigh(C)

        # sort from largest to smallest radius
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        R1_orig = np.sqrt(eigvals[0])
        R2_orig = np.sqrt(eigvals[1])

        # angle of the major axis
        angle_orig = np.arctan2(eigvecs[1, 0], eigvecs[0, 0])

        return xopt_orig, R1_orig, R2_orig, angle_orig, eigvecs

    # --------------------------------------------------------------------------

    def my_ellipse_err(self,lx):
        """
        Error function used to optimize ellipse axes and rotation.

        lx = [r1, r2, angle]
        """

        r1 = lx[0]
        r2 = lx[1]

        if r1 <= 0.0 or r2 <= 0.0:
            return np.inf

        angle = lx[2]

        err = 0.0
        cnt = 0.0

        for i in  range(0, 360, 5):
            t = i * np.pi / 180.0

            x = np.zeros(2)

            x[0] = r1 * np.cos(t) * np.cos(angle) - r2 * np.sin(t) * np.sin(angle) + self.xopt[0]
            x[1] = r1 * np.cos(t) * np.sin(angle) + r2 * np.sin(t) * np.cos(angle) + self.xopt[1]

            [ene, g, h] = self.eval_uv(x)

            if self.ts_ellipse == 1:
                z   = - self.thr * np.cos(2*t)
            else:
                z   = self.thr

            err += ((ene - self.xoptene) - z) ** 2
            cnt += 1.0

        return np.sqrt(err / cnt)

# =============================================================================
# detect basins
# =============================================================================

    def calculate_watershed(self):
        """
        Calculate watershed regions on the FES.

        The FES itself is treated as the topographic surface.
        Minima from points.txt are used as markers.
        """

        # Regions with energy greater or equal to zmax are treated as unsampled/high-energy boundary.
        self.basins_mask = np.isfinite(self.Z) & (self.Z < self.zmax)

        self.basins_markers = self.points_to_markers(self.basins_mask)

        self.basins_labels = watershed(
            self.Z,
            markers=self.basins_markers,
            mask=self.basins_mask,
            watershed_line=True,
        )

        return self.basins_labels, self.basins_markers, self.basins_mask

    # --------------------------------------------------------------------------

    def points_to_markers(self, mask):
        """
        Convert local minima coordinates to a marker image for watershed.

        Marker image convention:
            0 = no marker
            1, 2, ... = watershed seeds
        """

        markers = np.zeros(mask.shape, dtype=np.int32)

        used_pixels = {}

        for p in self.sp_optimized:

            if p['type'] != "S":
                continue

            marker_id = int(p['id'])
            ix = int(np.argmin(np.abs(self.x_grid - p['x'])))
            iy = int(np.argmin(np.abs(self.y_grid - p['y'])))

            if not mask[iy, ix]:
                print(
                    f"    WARNING: seed {p['id']} at ({p['x']:.4f}, {p['y']:.4f}) "
                    f"is outside the valid FES mask; searching nearest valid grid point."
                )

                valid_y, valid_x = np.where(mask)
                dist2 = (self.x_grid[valid_x] - p['x'])**2 + (self.y_grid[valid_y] - p['y'])**2
                nearest = np.argmin(dist2)

                iy = valid_y[nearest]
                ix = valid_x[nearest]

            if (iy, ix) in used_pixels:
                old_id = used_pixels[(iy, ix)]
                print(
                    f"    WARNING: seed {p['id']} maps to the same pixel as seed {old_id}; "
                    f"keeping the new marker at the same location."
                )

            markers[iy, ix] = marker_id
            used_pixels[(iy, ix)] = marker_id

        return markers

    # --------------------------------------------------------------------------

    def calculate_basins_ene_state(self):

        minima = [pt for pt in self.sp_optimized if pt["type"] == "S"]

        if len(minima) == 0:
            print("  WARNING: No minima found; basin state energies cannot be calculated.")
            return

        for pt in minima:
            region_id = int(pt['id'])
            dx = 1.0 / self.x_axis.npts
            dy = 1.0 / self.y_axis.npts
            mask = self.basins_labels == region_id
            weights = np.exp(-self.Z[mask] / (self.temp * rfac))
            Q =  np.sum(weights) * dx * dy
            ene_int = - self.temp * rfac * math.log(Q)
            pt['ene_state'] = ene_int

        # get minimum energy and calculate corrected energy value
        ene_min = min(pt['ene'] for pt in self.sp_optimized if pt['type'] == "S")
        ene_state_min = min(pt['ene_state'] for pt in self.sp_optimized if pt['type'] == "S")
        for pt in self.sp_optimized:
            if pt['type'] == "S":
                pt['ene0'] = pt['ene'] - ene_min
                pt['ene_state0'] = pt['ene_state'] - ene_state_min

# --------------------------------------------------------------------------

    def remove_transition_states_inside_minima_basins(self):
        """
        Remove transition-state stationary points that fall inside minima basins.

        Watershed labels are generated from local minima. A transition state should
        lie on a basin boundary. If a point of type ``T`` maps to a labelled basin
        pixel, it is inside the corresponding minimum basin and is removed from
        ``self.sp_optimized``. Points on watershed lines, outside the valid mask,
        or outside the grid are kept.

        Returns
        -------
        removed : list of dict
            Removed transition-state stationary points.
        """

        print("")
        print("# Detecting transition states inside minima basins ...")

        if self.basins_labels is None:
            raise RuntimeError("Basin labels are not available. Call calculate_watershed() first.")

        kept = []
        removed = []

        for item in self.sp_optimized:
            if item.get('type', "") != "T":
                kept.append(item)
                continue

            ix = int(np.argmin(np.abs(self.x_grid - item['x'])))
            iy = int(np.argmin(np.abs(self.y_grid - item['y'])))

            # If the nearest grid point is not the point itself because the point
            # is outside the plotting range, keep it rather than silently remove it.
            outside_grid = (
                item['x'] < self.x_axis.cvmin or item['x'] > self.x_axis.cvmax
                or item['y'] < self.y_axis.cvmin or item['y'] > self.y_axis.cvmax
            )

            if outside_grid:
                kept.append(item)
                continue

            basin_id = int(self.basins_labels[iy, ix])

            if basin_id > 0:
                removed.append(item)
                print(
                    "  Removing transition state inside minimum basin: "
                    f"label = {item['id']}, "
                    f"x = {item['x']:.6f}, "
                    f"y = {item['y']:.6f}, "
                    f"ene = {item['ene']:.6f}, "
                    f"basin = {basin_id}"
                )
            else:
                kept.append(item)

        self.sp_optimized = kept

        print(f"  Removed {len(removed)} transition state(s) inside minima basins.")
        print(f"  Remaining optimized stationary points: {len(self.sp_optimized)}")

        return removed

# ==============================================================================
# Command-line arguments
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

    parser = argparse.ArgumentParser(
        description="Detect and optimize stationary points on a 2D free-energy/potential-energy surface."
    )

    # --------------------------------------------------------------------------
    # CV1
    # --------------------------------------------------------------------------

    cv1group = parser.add_argument_group("The first collective variable (CV1) specification")

    cv1group.add_argument("--cv1label", type=str, default=r"cv1",
        help="Label of the first collective variable." )

    cv1group.add_argument("--cv1min", type=float, required=True,
        help="Minimum value of the first collective variable." )

    cv1group.add_argument("--cv1max", type=float, required=True,
        help="Maximum value of the first collective variable." )

    cv1group.add_argument("--cv1nbins", type=int, required=True,
        help="Number of grid bins/points for the first collective variable." )

    # --------------------------------------------------------------------------
    # CV2
    # --------------------------------------------------------------------------

    cv2group = parser.add_argument_group("The second collective variable (CV2) specification")

    cv2group.add_argument("--cv2label", type=str, default=r"cv2",
        help="Label of the second collective variable." )

    cv2group.add_argument("--cv2min", type=float, required=True,
        help="Minimum value of the second collective variable." )

    cv2group.add_argument("--cv2max", type=float, required=True,
        help="Maximum value of the second collective variable." )

    cv2group.add_argument("--cv2nbins", type=int, required=True,
        help="Number of grid bins/points for the second collective variable." )

    # --------------------------------------------------------------------------
    # Energy label
    # --------------------------------------------------------------------------

    enegroup = parser.add_argument_group("The energy axis specification")

    enegroup.add_argument("--enelabel", type=str, default=r"${\Delta}G [kcal/mol]$",
        help="Energy label." )

    enegroup.add_argument("--zmax", type=float, required=True,
        help="Maximum energy value considered." )

    enegroup.add_argument("--contour_spacing", type=float, default=1.0,
        help="Contour spacing." )

    # --------------------------------------------------------------------------
    # RBF interpolation
    # --------------------------------------------------------------------------

    rbfgroup = parser.add_argument_group("The RBF (Radial Basis Function) interpolation specification")

    rbfgroup.add_argument("--cv1nrbfs", type=int, default=20,
        help="Number of RBFs for the first collective variable." )

    rbfgroup.add_argument("--cv2nrbfs", type=int, default=20,
        help="Number of RBFs for the second collective variable." )

    rbfgroup.add_argument("--rbfwidthmode", type=str, default="static",
        help="RBF width mode: static, gridsearch, optimize." )

    rbfgroup.add_argument("--rbfsx", type=float, default=1.5,
        help="Width factor for CV1 in the RBF static width mode." )

    rbfgroup.add_argument("--rbfsy", type=float, default=1.5,
        help="Width factor for CV2 in the RBF static width mode." )

    rbfgroup.add_argument("--rcond", type=float, default=1.0e-9,
        help="SVD cutoff for RBF fitting." )

    # --------------------------------------------------------------------------
    # Files
    # --------------------------------------------------------------------------

    filegroup = parser.add_argument_group("The input/output files specification")

    filegroup.add_argument("--input-fes", type=str, dest="fname_input_fes", required=True,
        help="Input FES filename." )
    
    filegroup.add_argument("--input-fes-x-column", type=int, default=1,
        help="Index of x-column in the input FES file." )
    
    filegroup.add_argument("--input-fes-y-column", type=int, default=2,
        help="Index of y-column in the input FES file." )
    
    filegroup.add_argument("--input-fes-e-column", type=int, default=3,
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
    # Other setup
    # --------------------------------------------------------------------------

    sysgroup = parser.add_argument_group("The system specification")

    sysgroup.add_argument(
        "--actions",
        type=str,
        default="all",
        help="Comma separated lists of actions: \n"
             "* all - \n"
             "* loadfes \n"
             "* optrbfs \n"
             "* calcsurfs \n"
             "* loadguess \n"
             "* guess \n"
             "* loadopts \n"
             "* optimize \n"
             "* basins."
    )

    sysgroup.add_argument( "--temp", type=float, default=300.0,
        help="Thermodynamic temperature." )

    sysgroup.add_argument( "--thrfac", type=float, default=0.25,
        help="Energy threshold factor for defining ellipse representing a stationary point (thrfac * kB * T)." )

    sysgroup.add_argument( "--random_seed", type=int, default=None,
        help="Random generator seed." )

    # --------------------------------------------------------------------------
    # Treshold
    # --------------------------------------------------------------------------

    thresholdgroup = parser.add_argument_group("Treshold specification")
    
    thresholdgroup.add_argument( "--sng-cutoff", type=float, default=100.0,
        help="Cut-off value of SNG to consider a stationary point." )

    thresholdgroup.add_argument( "--sng-min-size", type=int, default=1,
        help="Minimum number of bins to consider a SNG region as a stationary point." )
    
    thresholdgroup.add_argument( "--trust-r0", type=float, default=0.01,
        help="Trust radius in scaled coordinates to locate stationary points." )
    
    thresholdgroup.add_argument( "--minslam", type=float, default=5.0,
        help="Minimal allowed value for scaled Hessian eigenvalue." )

    thresholdgroup.add_argument( "--maxslam", type=float, default=5000.0,
        help="Maximum allowed value for scaled Hessian eigenvalue." )

    thresholdgroup.add_argument( "--min-distance-uv", type=float, default=0.1,
        help="Minimum distance between two SPs to be considered as individual points." )

    thresholdgroup.add_argument( "--min-r", type=float, default=0.001,
        help="Minimum ellipse radius (scaled units)." )

    thresholdgroup.add_argument( "--max-sng", type=float, default=0.5,
        help="Maximum value of SNG for stationary points." )

    # --------------------------------------------------------------------------
    # Plots
    # --------------------------------------------------------------------------

    plotgroup = parser.add_argument_group("The graphical plot specification")

    plotgroup.add_argument( "--showrawfes", action="store_true",
        help="Show an interactive plot with the raw FES loaded (action: loadfes)." )

    plotgroup.add_argument( "--saverawfes", type=str, default="FigureFES-RAW.png",
        help="Save a plot with the raw FES loaded into a file (action: loadfes)." )

    plotgroup.add_argument( "--showrbffes", action="store_true",
        help="Show an interactive plot with the RBF interpolated FES (action: calcsurfs)." )

    plotgroup.add_argument( "--saverbffes", type=str, default="FigureFES-RBF.png",
        help="Save a plot with the RBF interpolated FES into a file (action: calcsurfs)." )

    plotgroup.add_argument( "--showsng", action="store_true",
        help="Show an interactive plot with the SNG regions (action: guess)." )

    plotgroup.add_argument( "--savesng", type=str, default="FigureFES-SNG.png",
        help="Save a plot with the SNG regions (action: guess)." )

    plotgroup.add_argument( "--showgpts", action="store_true",
        help="Show an interactive plot with the RBF interpolated FES and guessed stationary points (action: guess, loadguess)." )

    plotgroup.add_argument( "--savegpts", type=str, default="FigureFES-GPTS.png",
        help="Save a plot with the RBF interpolated FES and guessed stationary points into a file (action: guess, loadguess)." )

    plotgroup.add_argument( "--showopts", action="store_true",
        help="Show an interactive plot with the RBF interpolated FES and optimized stationary points (action: optimize, loadopts)." )

    plotgroup.add_argument( "--saveopts", type=str, default="FigureFES-OPTS.png",
        help="Save a plot with the RBF interpolated FES and optimized stationary points into a file (action: optimize, loadopts)." )

    plotgroup.add_argument( "--showbasins", action="store_true",
        help="Show an interactive plot with the RBF interpolated FES and minima basins (action: basins)." )

    plotgroup.add_argument( "--savebasins", type=str, default="FigureFES-Basins.png",
        help="Save a plot with the RBF interpolated FES and minima basins into a file (action: basins)." )

    plotgroup.add_argument('--figsize',type=parse_figsize,default=(6.4, 4.8),  # Default Matplotlib size fallback
        help="Figure size as 'width,height' in inches (default: 6.4,4.8)" )

    plotgroup.add_argument( "--dpi", type=int, default=300,
        help="Resolution for plot figures." )

    return parser.parse_args()

# ==============================================================================
# Action handlers
# ==============================================================================

def load_fes(args,surf):
    print("")
    print(f"# Load FES: {args.fname_input_fes}")
    surf.load(args.fname_input_fes,xcolumn=args.input_fes_x_column,ycolumn=args.input_fes_y_column,ecolumn=args.input_fes_e_column)

    if args.showrawfes == True or args.saverawfes is not None:
        surf.plot_raw(show=args.showrawfes,save=args.saverawfes,figsize=args.figsize,dpi=args.dpi)

# ------------------------------------------------------------------------------

def opt_rbfs(args,surf):
    print(f"")
    print(f"# Optimize RBF ...")
    print(f"  Width mode: {args.rbfwidthmode}")

    if args.rbfwidthmode == "static":
        print(f"  Sx:         {args.rbfsx:10.3f}")
        print(f"  Sy:         {args.rbfsy:10.3f}")
        surf.fit(sx=args.rbfsx,sy=args.rbfsy,rcond=args.rcond)
    elif args.rbfwidthmode == "gridsearch":
        surf.fit_width_grid(rcond=args.rcond)
    elif args.rbfwidthmode == "optimize":
        surf.fit_width_optimize(rcond=args.rcond)
    else:
        raise NameError(f"Unsupported action: {args.rbfwidthmode:s}")

    print(f"  RMSE:       {surf.fit_rmse:10.3f}")

# ------------------------------------------------------------------------------

def cal_surfs(args,surf):
    print("")
    print("# Calculate ENE and SNG surfaces ...")
    surf.calc_ene_and_sng()

    if args.showrbffes == True or args.saverbffes is not None:
        surf.plot_rbf(show=args.showrbffes,save=args.saverbffes,figsize=args.figsize,dpi=args.dpi)

# ------------------------------------------------------------------------------

def load_sp_guesses(args,surf):

    print("")
    print(f"# Load initial stationary point guesses (2nd column -> x, 3rd column -> y) ...")
    print(f"  File name: {args.fname_sp_guesses}")

    surf.sp_guesses = []

    with open(args.fname_sp_guesses, "r", encoding="utf-8") as fin:
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

            if len(fields) >= 2:
                indx = len(surf.sp_guesses) + 1
                surf.sp_guesses.append({
                    'id': int(indx),
                    'x': float(fields[1]),
                    'y': float(fields[2]),
                })
            else:
                raise ValueError("Stationary point file must contain at least 2 columns.")

    print(f"  Number of loaded stationary point guesses: {len(surf.sp_guesses)}")

    print("")
    print("#       ID          X          Y")
    print("# -------- ---------- ----------")

    for pt in surf.sp_guesses:
        print(f"{pt['id']:10d} {pt['x']:10.3f} {pt['y']:10.3f}")

    if args.showgpts == True or args.savegpts is not None:
        surf.plot_rbf_with_sp_guesses(show=args.showgpts,save=args.savegpts,figsize=args.figsize,dpi=args.dpi)

# ------------------------------------------------------------------------------

def sp_guess(args,surf):
    print("")
    print("# Detect stationary points ...")
    surf.detect_sp_regions_connected(
        eps=args.sng_cutoff,
        connectivity=1,
        representative="minabs",
        min_size=args.sng_min_size,
    )

    print(f"  Number of detected stationary points: {len(surf.sp_guesses)}")

    print("")
    print(f"# Found stationary points ...")
    print(f"  Saved as: {args.fname_sp_guesses}")
    with open(args.fname_sp_guesses, "w") as fout_spg:
        fout_spg.write("# Automatically detected stationary-point guesses\n")
        fout_spg.write("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)       Size\n")
        fout_spg.write("# -------- ---------- ---------- ---------- ---------- ----------\n")

        for pt in surf.sp_guesses:
            fout_spg.write(
                f"{pt['id']:10d} "
                f"{pt['x']:10.3f} "
                f"{pt['y']:10.3f} "
                f"{pt['ene']:10.3f} "
                f"{pt['sng']:10.3f} "
                f"{pt['size']:10d}\n"
            )

    print("")
    print("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)       Size")
    print("# -------- ---------- ---------- ---------- ---------- ----------")

    for pt in surf.sp_guesses:
        print(f"{pt['id']:10d} {pt['x']:10.3f} {pt['y']:10.3f} {pt['ene']:10.1f} {pt['sng']:10.1f} {pt["size"]:10d}")

    if args.showgpts == True or args.savegpts is not None:
        surf.plot_rbf_with_sp_guesses(show=args.showgpts,save=args.savegpts,figsize=args.figsize,dpi=args.dpi)

    if args.showsng == True or args.savesng is not None:
        surf.plot_uv_sng(show=args.showsng,save=args.savesng,figsize=args.figsize,dpi=args.dpi,zmax=args.sng_cutoff)

# ------------------------------------------------------------------------------

def load_sp_optimized(args,surf):

    print("")
    print(f"# Load optimized stationary points ...")
    print(f"  File name: {args.fname_sp_optimized}")

    surf.sp_optimized = []

    with open(args.fname_sp_optimized, "r", encoding="utf-8") as fin:
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

            if len(fields) >= 9:
                surf.sp_optimized.append({
                    'id': int(fields[0]),
                    'x': float(fields[1]),
                    'y': float(fields[2]),
                    'ene': float(fields[3]),
                    'sng': float(fields[4]),
                    'r1': float(fields[5]),
                    'r2': float(fields[6]),
                    'angle': float(fields[7]),
                    'type': str(fields[8])
                })
            else:
                raise ValueError("Stationary point file must contain at least 9 columns.")

    print(f"  Number of loaded optimized stationary points: {len(surf.sp_optimized)}")

    print("")
    print("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)         R1         R2      Angle       Type")
    print("# -------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        print(f"{opt['id']:10d} {opt['x']:10.3f} {opt['y']:10.3f} {opt['ene']:10.3f} {opt['sng']:10.1f} {opt['r1']:10.3f} {opt['r2']:10.3f} {opt['angle']:10.1f} {opt['type']:>10}")
    print("")

    if args.showopts == True or args.saveopts is not None:
        surf.plot_rbf_with_sp_optimized(show=args.showopts,save=args.saveopts,figsize=args.figsize,dpi=args.dpi)

# ------------------------------------------------------------------------------

def opt_sps(args,surf):
    print("")
    print("# Optimize stationary points ...")

    for gpt in surf.sp_guesses:
        print(f"")
        print(f"# ==============================================================================")
        print(f">>> Processing initial guess: x = {gpt['x']}, y = {gpt['y']}")
        surf.find_sp(gpt['x'], gpt['y'],trust_r0=args.trust_r0)

    print("# ==============================================================================")

    print(f"")
    print("# Cleaning stationary points ...")

    print("")
    print("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)         R1         R2      Angle Type        sL1        sL2        sR1        sR2")
    print("# -------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        print(f"{opt['id']:10d} {opt['x']:10.3f} {opt['y']:10.3f} {opt['ene']:10.3f} {opt['sng']:10.1f} {opt['r1']:10.3f} {opt['r2']:10.3f} {opt['angle']:10.1f} {opt['type']:>4} {opt['sl1']:10.1f} {opt['sl2']:10.1f} {opt['sr1']:10.3f} {opt['sr2']:10.3f}")

    surf.remove_sng_outliers(max_sng=args.max_sng)
    surf.remove_stationary_point_outliers(minslam=args.minslam,maxslam=args.maxslam)
    surf.remove_small_ellipses(minr=args.min_r)
    surf.remove_near_duplicate_sp_optimized(min_distance_uv=args.min_distance_uv)

    print(f"")
    print(f"# Optimized stationary points ...")
    print(f"  Saved as:                        {args.fname_sp_optimized}")
    print(f"  Number of all stationary points: {len(surf.sp_optimized)}")
    print(f"  Number of minima:                {sum(1 for item in surf.sp_optimized if item['type'] == "S")}")
    print(f"  Number of transition states:     {sum(1 for item in surf.sp_optimized if item['type'] == "T")}")
    print(f"  Number of maxima:                {sum(1 for item in surf.sp_optimized if item['type'] == "M")}")

    with open(args.fname_sp_optimized, "w") as fout_spo:
        fout_spo.write("# Optimized stationary points\n")
        fout_spo.write("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)         R1         R2      Angle       Type\n")
        fout_spo.write("# -------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n")

        for opt in surf.sp_optimized:
            fout_spo.write(
                f"{opt['id']:10d} "
                f"{opt['x']:10.3f} "
                f"{opt['y']:10.3f} "
                f"{opt['ene']:10.3f} "
                f"{opt['sng']:10.3f} "
                f"{opt['r1']:10.3f} "
                f"{opt['r2']:10.3f} "
                f"{opt['angle']:10.3f} "
                f"{opt['type']:>10s}\n"
            )

    print("")
    print("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)         R1         R2      Angle       Type")
    print("# -------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        print(f"{opt['id']:10d} {opt['x']:10.3f} {opt['y']:10.3f} {opt['ene']:10.3f} {opt['sng']:10.1f} {opt['r1']:10.3f} {opt['r2']:10.3f} {opt['angle']:10.1f} {opt['type']:>10}")

    if args.showopts == True or args.saveopts is not None:
        surf.plot_rbf_with_sp_optimized(show=args.showopts,save=args.saveopts,figsize=args.figsize,dpi=args.dpi)

# -------------------------------------------------------------------------

def find_basins(args,surf):
    print("")
    print("# Find minimum basins ...")

    surf.calculate_watershed()
    surf.remove_transition_states_inside_minima_basins()
    surf.calculate_basins_ene_state()

    print(f"  Number of minimum basins:        {sum(1 for item in surf.sp_optimized if item['type'] == "S")}")

    with open(args.fname_sp_basins, "w") as fout_spo:
        fout_spo.write("# Basin stationary points\n")
        fout_spo.write("#       ID          X          Y   ENE(X,Y)   SNG(X,Y)         R1         R2      Angle       Type\n")
        fout_spo.write("# -------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n")

        for opt in surf.sp_optimized:
            fout_spo.write(
                f"{opt['id']:10d} "
                f"{opt['x']:10.3f} "
                f"{opt['y']:10.3f} "
                f"{opt['ene']:10.3f} "
                f"{opt['sng']:10.3f} "
                f"{opt['r1']:10.3f} "
                f"{opt['r2']:10.3f} "
                f"{opt['angle']:10.3f} "
                f"{opt['type']:>10s}\n"
            )

    print("")
    print("#       ID Type          X          Y A=ENE(X,Y)    A0(X,Y)   A(state)  A0(state)")
    print("# -------- ---- ---------- ---------- ---------- ---------- ---------- ----------")
    for opt in surf.sp_optimized:
        if opt['type'] == "S":
            print(f"{opt['id']:10d} {opt['type']:>4} {opt['x']:10.3f} {opt['y']:10.3f} {opt['ene']:10.3f} {opt['ene0']:10.3f} {opt['ene_state']:10.3f} {opt['ene_state0']:10.3f}")
    print("")

    with open(args.fname_basins, "w") as fout_spo:
        fout_spo.write("# Basins\n")
        fout_spo.write("#       ID Type          X          Y A=ENE(X,Y)    A0(X,Y)   A(state)  A0(state)\n")
        fout_spo.write("# -------- ---- ---------- ---------- ---------- ---------- ---------- ----------\n")

        for opt in surf.sp_optimized:
            if opt['type'] == "S":
                fout_spo.write(
                    f"{opt['id']:10d} "
                    f"{opt['type']:>4s} "
                    f"{opt['x']:10.3f} "
                    f"{opt['y']:10.3f} "
                    f"{opt['ene']:10.3f} "
                    f"{opt['ene0']:10.3f} "
                    f"{opt['ene_state']:10.3f} "
                    f"{opt['ene_state0']:10.3f}\n"
                )

    if args.showbasins == True or args.savebasins is not None:
        surf.plot_rbf_with_watershed_basins(show=args.showbasins,save=args.savebasins,figsize=args.figsize,dpi=args.dpi)

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":

    print("")
    print("# ==============================================================================")
    print("#                        *** Analyze 2D Energy Surface ***                     #")
    print("#          The analyse-2D-surface utility is part of the PMFLib toolkit.       #")
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
    print("# Initializing 2D energy surface ...")

    surf = EnergySurface2D(
        cv1min=args.cv1min,
        cv1max=args.cv1max,
        cv1_nx=args.cv1nrbfs,
        cv1_px=args.cv1nbins,
        cv1_label=args.cv1label,

        cv2min=args.cv2min,
        cv2max=args.cv2max,
        cv2_nx=args.cv2nrbfs,
        cv2_px=args.cv2nbins,
        cv2_label=args.cv2label,

        zmax=args.zmax,
        ene_label=args.enelabel,
        contour_spacing=args.contour_spacing,
        temp=args.temp,
        thrfac=args.thrfac,

        random_seed=args.random_seed,
    )

    actions = [item.strip() for item in args.actions.split(",")]

    if "all" in actions:
        actions = ["loadfes", "optrbfs", "calcsurfs", "guess", "optimize", "basins" ]

    print(f"  Actions: {", ".join(actions)}")

    for action in actions:
        if action == "loadfes":
            load_fes(args,surf)
        elif action == "optrbfs":
            opt_rbfs(args,surf)
        elif action == "calcsurfs":
            cal_surfs(args,surf)
        elif action == "loadguess":
            load_sp_guesses(args,surf)
        elif action == "guess":
            sp_guess(args,surf)
        elif action == "loadopts":
            load_sp_optimized(args,surf)
        elif action == "optimize":
            opt_sps(args,surf)
        elif action == "basins":
            find_basins(args,surf)
        else:
            raise NameError(f"Unsupported action: {action:s}")

    print("")

# ==============================================================================



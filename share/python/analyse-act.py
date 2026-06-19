#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================

def read_columns(filename, columns):
    """
    Read selected columns from a text file.

    Lines starting with '#' are ignored.
    Column indices are zero-based internally.
    """

    data = []

    with open(filename, "r") as fin:
        for line in fin:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split()

            try:
                row = [float(fields[col]) for col in columns]
            except IndexError:
                raise ValueError(
                    f"Line has fewer columns than requested: '{line}'"
                )
            except ValueError:
                raise ValueError(
                    f"Cannot convert line to numbers: '{line}'"
                )

            data.append(row)

    if len(data) == 0:
        raise ValueError("No numerical data were read from the file.")

    return np.asarray(data)

# ==============================================================================

def autocorrelation_fft(x):
    """
    Calculate normalized autocorrelation function using FFT.

    The returned ACF satisfies acf[0] = 1.
    """

    x = np.asarray(x, dtype=float)
    n = len(x)

    if n < 2:
        raise ValueError("At least two data points are required.")

    x = x - np.mean(x)

    variance = np.var(x)

    if variance == 0.0:
        raise ValueError("The variance of the time series is zero.")

    # Zero padding to avoid circular correlation
    nfft = 1 << (2 * n - 1).bit_length()

    fx = np.fft.fft(x, n=nfft)
    acf = np.fft.ifft(fx * np.conjugate(fx)).real[:n]

    # Unbiased normalization by the number of overlapping points
    acf /= np.arange(n, 0, -1)

    # Normalize so that ACF(0) = 1
    acf /= acf[0]

    return acf

# ==============================================================================

def integrated_autocorrelation_time(acf, cutoff="first-negative"):
    """
    Calculate integrated autocorrelation time from the ACF.

    In units of the sampling interval:

        tau_int = 1 + 2 * sum_{k=1}^{K} ACF(k)

    The cutoff K can be determined by the first negative ACF value.
    """

    if cutoff != "first-negative":
        raise ValueError("Only cutoff='first-negative' is currently supported.")

    tau = 1.0

    for k in range(1, len(acf)):
        if acf[k] <= 0.0:
            break
        tau += 2.0 * acf[k]

    return tau

# ==============================================================================

def save_acf(filename, acfs, dt_fs):
    """
    Save ACF values for all analysed columns.

    Output columns:
        index, time_fs, acf_col1, acf_col2, ...
    """

    n = len(acfs[0])
    time = np.arange(n) * dt_fs

    table = [np.arange(n), time]

    for acf in acfs:
        table.append(acf)

    table = np.column_stack(table)

    header = "index time_fs " + " ".join(
        f"acf_col{i+1}" for i in range(len(acfs))
    )

    np.savetxt(filename, table, header=header)

# ==============================================================================

def plot_acf(acfs, columns, dt_fs, time_unit, show, filename, figsize, dpi):
    """
    Plot ACFs for all analysed columns.
    """

    plt.figure(figsize=figsize)

    if time_unit == "fs":
        x = np.arange(len(acfs[0])) * dt_fs
        xlabel = "Time / fs"
    else:
        x = np.arange(len(acfs[0]))
        xlabel = "Sampling step"

    for acf, col in zip(acfs, columns):
        plt.plot(x, acf, label=f"column {col + 1}")

    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel(xlabel)
    plt.ylabel("ACF")
    plt.legend()
    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, dpi=dpi)

    if show:
        plt.show()

    plt.close()

# ==============================================================================

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

# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate autocorrelation time for selected columns."
    )

    parser.add_argument("--input",type=str,required=True,
        help="Input text file with numerical data." )

    parser.add_argument("-c", "--columns",required=True,
        help="Columns to analyse, using one-based indices. Example: -c 1,2,5"
    )

    parser.add_argument("--dt",type=float,default=1.0,
        help="Sampling period in steps or femtoseconds." )

    parser.add_argument("--time-unit", choices=["steps", "fs"],default="steps",
        help="Unit used for reporting autocorrelation time. Use 'steps' for sampling-time units or 'fs' for femtoseconds.")

    parser.add_argument("--save-acf",default=None,
        help="Output file for saving the ACF." )
    
    parser.add_argument("--save-acf-plot",default=None,
        help="Output file for saving the ACF." )

    parser.add_argument("--show-acf-plot",action="store_true",
        help="Show ACF plot." )
    
    parser.add_argument('--figsize',type=parse_figsize,default=(6.4, 4.8),  # Default Matplotlib size fallback
        help="Figure size as 'width,height' in inches (default: 6.4,4.8)" )
    
    parser.add_argument("--dpi", type=int, default=300,
        help="Resolution for plot figures.")

    args = parser.parse_args()

    # Convert one-based column indices to zero-based indices
    columns = [int(c.strip()) - 1 for c in args.columns.split(",")]

    if any(c < 0 for c in columns):
        raise ValueError("Column indices must be one-based and positive.")

    data = read_columns(args.input, columns)

    acfs = []
    taus_steps = []

    for i in range(data.shape[1]):
        acf = autocorrelation_fft(data[:, i])
        tau_steps = integrated_autocorrelation_time(acf)

        acfs.append(acf)
        taus_steps.append(tau_steps)

    for col, tau_steps in zip(columns, taus_steps):
        if args.time_unit == "fs":
            tau = tau_steps * args.dt
            print(f"column#{col + 1:03d}: tau_int = {tau:15.1f} fs")
        else:
            print(f"column#{col + 1:03d}: tau_int = {tau_steps:15.1f} steps")

    if args.save_acf is not None:
        save_acf(args.save_acf, acfs, args.dt, )

    if args.save_acf_plot is not None or args.show_acf_plot == True:
        plot_acf(acfs, columns, args.dt, args.time_unit, args.show_acf_plot, args.save_acf_plot, args.figsize, args.dpi)


if __name__ == "__main__":
    main()



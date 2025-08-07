"""
Plotting utilities for EVPA rotation analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from .bayesian_blocks import compute_bayesian_blocks, compute_local_extrema


# Configure plotting style
def setup_plot_style():
    """Configure matplotlib plotting style."""
    try:
        plt.style.use(['science', 'notebook', 'grid'])
    except OSError:
        # Fallback if scienceplots is not available
        plt.style.use('default')
    
    plt.rcParams.update({
        "axes.grid": True,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"]
    })
    
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['figure.dpi'] = 100


def plot_bayesian_blocks_with_extrema(data, source, p0=0.001, show_uncertainty=True):
    """
    Plots adjusted EVPA data with Bayesian Blocks and local extrema.
    Fills between extrema with blue for decreasing and orange for increasing EVPA.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data containing source observations
    source : str
        Source name to plot
    p0 : float, default=0.001
        Prior probability for Bayesian Blocks
    show_uncertainty : bool, default=True
        Whether to show uncertainty regions
    """
    from .loader import prepare_source_data
    
    df_source = prepare_source_data(data, source)
    
    for seg, group in df_source.groupby('segment'):
        if len(group) <= 4:
            continue
            
        edges, bin_means, bin_errors = compute_bayesian_blocks(group, p0=p0)

        # Skip segment if all means are NaN (no valid bins)
        if np.all(np.isnan(bin_means)):
            continue

        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        maxima_dict, minima_dict = compute_local_extrema(bin_means, bin_centers)

        # Combine extrema for filling
        extrema = sorted(
            [(mjd, evpa, 'max') for mjd, evpa in zip(maxima_dict['MJD'], maxima_dict['EVPA'])] +
            [(mjd, evpa, 'min') for mjd, evpa in zip(minima_dict['MJD'], minima_dict['EVPA'])],
            key=lambda x: x[0]
        )

        y_values = np.concatenate([np.array(bin_means), [bin_means[-1]]])

        plt.figure(figsize=(10, 6))
        
        # Plot data points
        plt.errorbar(group['MJD'], group['adjusted_EVPA'],
                     yerr=group['err_EVPA[deg]'], fmt='o', capsize=3,
                     color='black', ecolor='#444444', label="Data Points")

        # Plot Bayesian Blocks
        plt.step(edges, y_values, where='post', color='red', 
                linewidth=2, label='Bayesian Blocks')
        
        # Plot extrema
        plt.scatter(maxima_dict['MJD'], maxima_dict['EVPA'], 
                   color='blue', marker='^', s=100, label='Local Maxima')
        plt.scatter(minima_dict['MJD'], minima_dict['EVPA'], 
                   color='green', marker='v', s=100, label='Local Minima')

        # Vertical lines at extrema
        for mjd, evpa, kind in extrema:
            plt.axvline(x=mjd, color='blue' if kind == 'max' else 'green', 
                       linestyle='--', alpha=0.6)

        y_min, y_max = plt.ylim()

        # Fill between extrema (increasing/decreasing regions)
        for i in range(len(extrema) - 1):
            mjd1, evpa1, _ = extrema[i]
            mjd2, evpa2, _ = extrema[i + 1]
            fill_color = 'orange' if evpa2 > evpa1 else 'skyblue'
            plt.fill_betweenx([y_min, y_max], mjd1, mjd2, 
                             color=fill_color, alpha=0.1)

        # Uncertainty region
        if show_uncertainty:
            for i in range(len(bin_means)):
                if not np.isnan(bin_means[i]):
                    plt.fill_between([edges[i], edges[i+1]],
                                   [bin_means[i] - bin_errors[i]] * 2,
                                   [bin_means[i] + bin_errors[i]] * 2,
                                   step='mid', color='lightcoral', alpha=0.3)

        plt.ylim(y_min, y_max)
        plt.title(f"{source} - Segment {seg} (Bayesian Blocks on Adjusted EVPA)\n"
                 f"(Consecutive points with gap ≤ 30 days, > 4 data points)")
        plt.xlabel("Modified Julian Date (MJD)")
        plt.ylabel("Adjusted EVPA [deg]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


def plot_rotation_event(sub_data, period_data, source_name, period_index, 
                       start_edge, end_edge, amplitude, t_pval, binom_p,
                       edges=None, bin_means=None, bin_errors=None,
                       output_dir=None):
    """
    Plot a rotation event with context and statistics.
    
    Parameters:
    -----------
    sub_data : pd.DataFrame
        Data for the rotation interval
    period_data : pd.DataFrame
        Full period data for context
    source_name : str
        Name of the source
    period_index : int
        Period identifier
    start_edge, end_edge : float
        Time boundaries of rotation
    amplitude : float
        EVPA amplitude of rotation
    t_pval : float
        t-test p-value
    binom_p : float
        Binomial test p-value
    edges, bin_means, bin_errors : arrays, optional
        Bayesian blocks data for full period
    output_dir : str, optional
        Directory to save plot
    """
    plt.figure(figsize=(10, 6))
    
    # Plot all period data (faint)
    if period_data is not None:
        plt.errorbar(period_data['MJD'], period_data['adjusted_EVPA'],
                     yerr=period_data['err_EVPA[deg]'], fmt='o',
                     capsize=3, color='k', alpha=0.3, label='Period data')
    
    # Plot rotation data
    plt.errorbar(sub_data['MJD'], sub_data['adjusted_EVPA'],
                 yerr=sub_data['err_EVPA[deg]'], fmt='o',
                 capsize=3, color='k', ecolor='k', label='Rotation data')

    # Plot full period Bayesian blocks if provided
    if edges is not None and bin_means is not None:
        full_y = np.concatenate([bin_means, [bin_means[-1]]])
        plt.step(edges, full_y, where='post',
                 color='#333333', alpha=0.7, linewidth=1.5, 
                 label='Bayesian blocks (full)')
        
        if bin_errors is not None:
            lower = np.concatenate([bin_means - bin_errors,
                                  [bin_means[-1] - bin_errors[-1]]])
            upper = np.concatenate([bin_means + bin_errors,
                                  [bin_means[-1] + bin_errors[-1]]])
            plt.fill_between(edges, lower, upper,
                           step='post', alpha=0.2, color='#333333')

    # Set limits
    plt.xlim(start_edge, end_edge)
    y_margin = amplitude * 0.2
    plt.ylim(sub_data['adjusted_EVPA'].min() - y_margin,
             sub_data['adjusted_EVPA'].max() + y_margin)

    plt.xlabel('Modified Julian Date (days)', fontsize=14)
    plt.ylabel('Corrected EVPA (°)', fontsize=14)
    plt.title(
        f'Rotation Event for {source_name} | Period {period_index} | '
        f'MJD {start_edge:.1f}–{end_edge:.1f} | ΔEVPA = {amplitude:.1f}°',
        fontsize=16
    )
    
    # Add statistics text box
    plt.text(
        0.05, 0.95,
        f"$p_{{\\mathrm{{t-test}}}} = {t_pval:.2e}$\n"
        f"$p_{{\\mathrm{{binomial}}}} = {binom_p:.2e}$",
        transform=plt.gca().transAxes,
        fontsize=12,
        ha='left', va='top',
        bbox=dict(
            facecolor='white',
            edgecolor='black',
            alpha=0.8,
            pad=5
        )
    )
    
    plt.legend(loc='best', frameon=True, framealpha=0.8)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    
    # Save if output directory provided
    if output_dir:
        import os
        filename = f"{source_name.replace(' ', '_')}_p{period_index}_{start_edge:.1f}_{end_edge:.1f}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_source_overview(data, source_name):
    """
    Create an overview plot showing all segments for a source.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Full dataset
    source_name : str
        Name of the source to plot
    """
    from .loader import prepare_source_data
    
    df_source = prepare_source_data(data, source_name)
    
    fig, axes = plt.subplots(len(df_source['segment'].unique()), 1, 
                            figsize=(12, 4 * len(df_source['segment'].unique())),
                            sharex=True)
    
    if len(df_source['segment'].unique()) == 1:
        axes = [axes]
    
    for i, (seg, group) in enumerate(df_source.groupby('segment')):
        ax = axes[i]
        
        ax.errorbar(group['MJD'], group['adjusted_EVPA'],
                   yerr=group['err_EVPA[deg]'], fmt='o-',
                   capsize=3, color='black', alpha=0.7)
        
        ax.set_title(f"Segment {seg} ({len(group)} points)")
        ax.set_ylabel("Adjusted EVPA [deg]")
        ax.grid(True)
    
    axes[-1].set_xlabel("Modified Julian Date (MJD)")
    plt.suptitle(f"EVPA Overview for {source_name}", fontsize=16)
    plt.tight_layout()
    plt.show()


# Initialize plotting style when module is imported
setup_plot_style()

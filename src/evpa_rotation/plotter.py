"""
Automatic plotting utilities for EVPA rotation analysis.
Plots are saved to disk and figures are closed automatically to allow batch processing.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from .bayesian_blocks import compute_bayesian_blocks, compute_local_extrema

# Configure plotting style
def setup_plot_style():
    try:
        plt.style.use(['science', 'notebook', 'grid'])
    except OSError:
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

setup_plot_style()


def save_and_close_plot(fig, filename, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_bayesian_blocks_with_extrema(data, source, output_dir="plots", p0=0.001, show_uncertainty=True):
    from .loader import prepare_source_data
    df_source = prepare_source_data(data, source)
    
    for seg, group in df_source.groupby('segment'):
        if len(group) <= 4:
            continue
        
        edges, bin_means, bin_errors = compute_bayesian_blocks(group, p0=p0)
        if np.all(np.isnan(bin_means)):
            continue
        
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        maxima_dict, minima_dict = compute_local_extrema(bin_means, bin_centers)

        extrema = sorted(
            [(mjd, evpa, 'max') for mjd, evpa in zip(maxima_dict['MJD'], maxima_dict['EVPA'])] +
            [(mjd, evpa, 'min') for mjd, evpa in zip(minima_dict['MJD'], minima_dict['EVPA'])],
            key=lambda x: x[0]
        )

        y_values = np.concatenate([np.array(bin_means), [bin_means[-1]]])
        fig = plt.figure(figsize=(10, 6))
        
        plt.errorbar(group['MJD'], group['adjusted_EVPA'],
                     yerr=group['err_EVPA[deg]'], fmt='o', capsize=3,
                     color='black', ecolor='#444444', label="Data Points")
        plt.step(edges, y_values, where='post', color='red', linewidth=2, label='Bayesian Blocks')
        plt.scatter(maxima_dict['MJD'], maxima_dict['EVPA'], color='blue', marker='^', s=100, label='Local Maxima')
        plt.scatter(minima_dict['MJD'], minima_dict['EVPA'], color='green', marker='v', s=100, label='Local Minima')
        for mjd, evpa, kind in extrema:
            plt.axvline(x=mjd, color='blue' if kind == 'max' else 'green', linestyle='--', alpha=0.6)

        y_min, y_max = plt.ylim()
        for i in range(len(extrema) - 1):
            mjd1, evpa1, _ = extrema[i]
            mjd2, evpa2, _ = extrema[i + 1]
            fill_color = 'orange' if evpa2 > evpa1 else 'skyblue'
            plt.fill_betweenx([y_min, y_max], mjd1, mjd2, color=fill_color, alpha=0.1)

        if show_uncertainty:
            for i in range(len(bin_means)):
                if not np.isnan(bin_means[i]):
                    plt.fill_between([edges[i], edges[i+1]],
                                     [bin_means[i] - bin_errors[i]] * 2,
                                     [bin_means[i] + bin_errors[i]] * 2,
                                     step='mid', color='lightcoral', alpha=0.3)

        plt.ylim(y_min, y_max)
        plt.title(f"{source} - Segment {seg} (Bayesian Blocks on Adjusted EVPA)")
        plt.xlabel("Modified Julian Date (MJD)")
        plt.ylabel("Adjusted EVPA [deg]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        filename = f"{source}_segment{seg}.png"
        save_and_close_plot(fig, filename, output_dir)


def plot_rotation_event(sub_data, period_data, source_name, period_index, 
                       start_edge, end_edge, amplitude, t_pval, binom_p,
                       edges=None, bin_means=None, bin_errors=None,
                       output_dir="plots"):
    fig = plt.figure(figsize=(10, 6))
    
    if period_data is not None:
        plt.errorbar(period_data['MJD'], period_data['adjusted_EVPA'],
                     yerr=period_data['err_EVPA[deg]'], fmt='o',
                     capsize=3, color='k', alpha=0.3, label='Period data')
    
    plt.errorbar(sub_data['MJD'], sub_data['adjusted_EVPA'],
                 yerr=sub_data['err_EVPA[deg]'], fmt='o',
                 capsize=3, color='k', ecolor='k', label='Rotation data')

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

    plt.xlim(start_edge, end_edge)
    y_margin = amplitude * 0.2
    plt.ylim(sub_data['adjusted_EVPA'].min() - y_margin,
             sub_data['adjusted_EVPA'].max() + y_margin)

    plt.xlabel('Modified Julian Date (days)', fontsize=14)
    plt.ylabel('Corrected EVPA (°)', fontsize=14)
    plt.title(f'Rotation Event for {source_name} | Period {period_index} | '
              f'MJD {start_edge:.1f}–{end_edge:.1f} | ΔEVPA = {amplitude:.1f}°', fontsize=16)
    
    plt.text(0.05, 0.95,
             f"$p_{{\\mathrm{{t-test}}}} = {t_pval:.2e}$\n"
             f"$p_{{\\mathrm{{binomial}}}} = {binom_p:.2e}$",
             transform=plt.gca().transAxes,
             fontsize=12, ha='left', va='top',
             bbox=dict(facecolor='white', edgecolor='black', alpha=0.8, pad=5))
    
    plt.legend(loc='best', frameon=True, framealpha=0.8)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    
    filename = f"{source_name.replace(' ', '_')}_rotation_p{period_index}_{start_edge:.1f}_{end_edge:.1f}.png"
    save_and_close_plot(fig, filename, output_dir)


def plot_source_overview(data, source_name, output_dir="plots"):
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
    
    filename = f"{source_name.replace(' ', '_')}_overview.png"
    save_and_close_plot(fig, filename, output_dir)

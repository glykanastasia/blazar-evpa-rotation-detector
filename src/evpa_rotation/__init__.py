"""
EVPA Rotation Analysis Package

A Python package for detecting and analyzing Electric Vector Position Angle (EVPA) 
rotations in polarimetric monitoring data of astronomical sources.

This package provides tools for:
- Loading and preprocessing polarimetric monitoring data
- Adjusting EVPA angles to resolve 180° ambiguity
- Performing Bayesian Blocks analysis on time series
- Detecting rotation events through statistical tests
- Visualizing results with publication-quality plots

Main Classes:
- RotationAnalyzer: Main analysis class for detecting rotation events

Main Functions:
- load_data: Load and preprocess CSV monitoring data
- adjust_angles: Resolve 180° EVPA ambiguity
- compute_bayesian_blocks: Perform Bayesian Blocks segmentation
- plot_bayesian_blocks_with_extrema: Create overview plots

Example usage:
    from evpa_rotation import RotationAnalyzer
    
    # Initialize analyzer
    analyzer = RotationAnalyzer(data_file="monitoring_data.csv")
    
    # Analyze all sources
    results = analyzer.analyze_all_sources()
    
    # Get summary statistics
    summary = analyzer.get_summary_statistics()
"""

from .main import RotationAnalyzer
from .loader import load_data, prepare_source_data, segment_data_by_gaps
from .angles import adjust_angles, compute_angle_statistics, circular_mean
from .bayesian_blocks import (
    compute_bayesian_blocks, 
    compute_local_extrema,
    analyze_extrema_pairs
)
from .plotter import (
    plot_bayesian_blocks_with_extrema,
    plot_rotation_event,
    plot_source_overview,
    setup_plot_style
)

__version__ = "1.0.0"
__author__ = "Anastasia Glykopoulou"
__email__ = "glykopoulouanastasia@gmail.com"

__all__ = [
    # Main classes
    "RotationAnalyzer",
    
    # Data loading and preprocessing
    "load_data",
    "prepare_source_data", 
    "segment_data_by_gaps",
    
    # Angle processing
    "adjust_angles",
    "compute_angle_statistics",
    "circular_mean",
    
    # Bayesian blocks analysis
    "compute_bayesian_blocks",
    "compute_local_extrema", 
    "analyze_extrema_pairs",
    
    # Plotting functions
    "plot_bayesian_blocks_with_extrema",
    "plot_rotation_event",
    "plot_source_overview",
    "setup_plot_style",
]

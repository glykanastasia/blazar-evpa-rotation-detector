"""
Main analysis module for EVPA rotation detection (corrected with error propagation).
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

from .loader import load_data, prepare_source_data
from .bayesian_blocks import (
    compute_bayesian_blocks, 
    compute_local_extrema,
    get_bin_boundaries_for_interval
)
from .plotter import plot_rotation_event

class RotationAnalyzer:
    """
    Main class for analyzing EVPA rotations in polarimetric monitoring data.
    """
    
    def __init__(self, data_file=None, data=None, output_dir="rotation_analysis"):
        if data is not None:
            self.data = data
        elif data_file is not None:
            self.data = load_data(data_file)
        else:
            raise ValueError("Either data_file or data must be provided")
        
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)
        
        self.rotation_data = []
        self.recorded_rotations = set()
        
    def analyze_source(self, source_name, p0=0.001, diff_threshold=0, 
                       t_test_threshold=0.05, binom_threshold=0.0625,
                       save_plots=True, show_plots=False):
        print(f"\nAnalyzing source: {source_name}")
        
        df_source = prepare_source_data(self.data, source_name)
        total_periods = df_source['segment'].nunique()
        print(f"Found {total_periods} observation periods")
        
        source_rotations = []
        
        for period_index in df_source['segment'].unique():
            period = df_source[df_source['segment'] == period_index].copy()
            if len(period) < 5:
                continue
            print(f"  Period {period_index}: {len(period)} points")
            
            period_rotations = self._analyze_period(
                period, source_name, period_index, p0, diff_threshold,
                t_test_threshold, binom_threshold, save_plots, show_plots
            )
            source_rotations.extend(period_rotations)
        
        return source_rotations
    
    def _analyze_period(self, period, source_name, period_index, p0, diff_threshold,
                        t_test_threshold, binom_threshold, save_plots, show_plots):
        
        period_rotations = []
        period = period.drop_duplicates(subset='MJD').sort_values('MJD')
        
        try:
            edges, bin_means, bin_errors = compute_bayesian_blocks(period, p0=p0)
        except Exception as e:
            print(f"    Bayesian blocks failed: {e}")
            return period_rotations
        
        bin_means = np.array(bin_means)
        bin_errors = np.array(bin_errors)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        
        maxima_dict, minima_dict = compute_local_extrema(bin_means, bin_centers)
        all_extrema = [(m, v, 'max') for m, v in zip(maxima_dict['MJD'], maxima_dict['EVPA'])]
        all_extrema += [(m, v, 'min') for m, v in zip(minima_dict['MJD'], minima_dict['EVPA'])]
        all_extrema.sort(key=lambda x: x[0])
        if len(all_extrema) < 2:
            return period_rotations
        
        existing_intervals = [
            (r['t_start'], r['t_end'])
            for r in self.rotation_data
            if r['Source'] == source_name
        ]
        
        for i in range(len(all_extrema) - 1):
            t_start, evpa_start, type_start = all_extrema[i]
            t_end, evpa_end, type_end = all_extrema[i + 1]
            
            if not ((type_start == 'max' and type_end == 'min') or
                    (type_start == 'min' and type_end == 'max')):
                continue
            
            if any(not (t_end <= s or t_start >= e) for s, e in existing_intervals):
                continue
            
            diff_evpa = abs(evpa_end - evpa_start)
            if diff_evpa < diff_threshold:
                continue
            
            start_edge, end_edge = get_bin_boundaries_for_interval(edges, t_start, t_end)
            sub_data = period[(period['MJD'] >= start_edge) & (period['MJD'] <= end_edge)].sort_values('MJD')
            if len(sub_data) < 4:
                continue
            
            try:
                edges_sub, means_sub, errs_sub = compute_bayesian_blocks(sub_data, p0=p0)
            except Exception as e:
                continue
            if (len(edges_sub) - 1) <= 2:
                continue
            
            # Statistical tests
            t_stat, t_pval = stats.ttest_1samp(sub_data['adjusted_EVPA'], popmean=180)
            n_above = (sub_data['adjusted_EVPA'] > 180).sum()
            binom_p = stats.binomtest(n_above, len(sub_data), p=0.5).pvalue
            
            # Amplitude and error
            amp = sub_data['adjusted_EVPA'].max() - sub_data['adjusted_EVPA'].min()
            emax = sub_data.loc[sub_data['adjusted_EVPA'].idxmax(), 'err_EVPA[deg]']
            emin = sub_data.loc[sub_data['adjusted_EVPA'].idxmin(), 'err_EVPA[deg]']
            amp_err = np.hypot(emax, emin)
            
            if amp < diff_threshold:
                continue
            
            rotation_key = (source_name, round(start_edge, 1), round(end_edge, 1))
            
            if (t_pval <= t_test_threshold and 
                binom_p <= binom_threshold and 
                amp >= diff_threshold and 
                rotation_key not in self.recorded_rotations):
                
                # ΔMJD error from start/end block uncertainties
                delta_mjd = end_edge - start_edge
                delta_mjd_err = np.hypot(bin_errors[0], bin_errors[-1])
                
                # Rotation rate
                rate = amp / delta_mjd if delta_mjd != 0 else 0
                rate_err = rate * np.hypot(amp_err/amp, delta_mjd_err/delta_mjd) if delta_mjd != 0 else 0
                
                rotation_event = {
                    'Source': source_name,
                    'AMPLITUDE': amp,
                    'AMPLITUDE_error': amp_err,
                    't_start': start_edge,
                    't_end': end_edge,
                    'ΔMJD': delta_mjd,
                    'ΔMJD_error': delta_mjd_err,
                    't-test p': t_pval,
                    'Binom test p': binom_p,
                    'Rate': rate,
                    'Rate_error': rate_err,
                    'n_points': len(sub_data),
                    'extrema_type': f"{type_start}→{type_end}"
                }
                
                self.rotation_data.append(rotation_event)
                self.recorded_rotations.add(rotation_key)
                existing_intervals.append((start_edge, end_edge))
                period_rotations.append(rotation_event)
                
                print(f"      ROTATION DETECTED: MJD {start_edge:.1f}-{end_edge:.1f}, "
                      f"Amp={amp:.1f}±{amp_err:.1f}°, ΔMJD={delta_mjd:.2f}±{delta_mjd_err:.2f}, "
                      f"t-p={t_pval:.2e}, binom-p={binom_p:.2e}")
                
                if save_plots or show_plots:
                    plot_output_dir = os.path.join(self.output_dir, "plots") if save_plots else None
                    plot_rotation_event(
                        sub_data, period, source_name, period_index,
                        start_edge, end_edge, amp, t_pval, binom_p,
                        edges, bin_means, bin_errors, plot_output_dir
                    )
                    if not show_plots:
                        import matplotlib.pyplot as plt
                        plt.close()
        
        return period_rotations
    
    def analyze_all_sources(self, sources=None, **kwargs):
        if sources is None:
            sources = self.data['J2000_name'].unique()
        
        all_rotations = []
        for source in tqdm(sources, desc="Analyzing sources"):
            try:
                source_rotations = self.analyze_source(source, **kwargs)
                all_rotations.extend(source_rotations)
            except Exception as e:
                print(f"Error analyzing {source}: {e}")
                continue
        
        if self.rotation_data:
            df = pd.DataFrame(self.rotation_data)
            # Format CSV with desired columns
            df_out = pd.DataFrame({
                'Source': df['Source'],
                'Amplitude (deg)': df.apply(lambda x: f"{x['AMPLITUDE']:.1f} ± {x['AMPLITUDE_error']:.1f}", axis=1),
                'Start (MJD)': df['t_start'],
                'End (MJD)': df['t_end'],
                'ΔMJD (days)': df.apply(lambda x: f"{x['ΔMJD']:.2f} ± {x['ΔMJD_error']:.2f}", axis=1),
                'p-t-test': df['t-test p'],
                'p-binom': df['Binom test p'],
                'Rate (deg/day)': df.apply(lambda x: f"{x['Rate']:.2f} ± {x['Rate_error']:.2f}", axis=1)
            })
            output_file = os.path.join(self.output_dir, "rotation_events.csv")
            df_out.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
            return df_out
        else:
            print("No rotation events detected.")
            return pd.DataFrame()
    
    def get_summary_statistics(self):
        if not self.rotation_data:
            return {"message": "No rotation events detected"}
        
        df = pd.DataFrame(self.rotation_data)
        summary = {
            "total_rotations": len(df),
            "unique_sources": df['Source'].nunique(),
            "sources_with_rotations": df['Source'].unique().tolist(),
            "amplitude_stats": {
                "mean": df['AMPLITUDE'].mean(),
                "std": df['AMPLITUDE'].std(),
                "min": df['AMPLITUDE'].min(),
                "max": df['AMPLITUDE'].max(),
                "median": df['AMPLITUDE'].median()
            },
            "duration_stats": {
                "mean": df['ΔMJD'].mean(),
                "std": df['ΔMJD'].std(),
                "min": df['ΔMJD'].min(),
                "max": df['ΔMJD'].max(),
                "median": df['ΔMJD'].median()
            },
            "rotations_per_source": df['Source'].value_counts().to_dict()
        }
        return summary

def main():
    import argparse
    parser = argparse.ArgumentParser(description="EVPA Rotation Analysis")
    parser.add_argument("--data", required=True, help="Path to CSV data file")
    parser.add_argument("--output", default="rotation_analysis", help="Output directory")
    parser.add_argument("--sources", nargs="+", help="Specific sources to analyze (default: all)")
    parser.add_argument("--p0", type=float, default=0.001, help="Prior probability for Bayesian blocks")
    parser.add_argument("--min-amp", type=float, default=0, help="Minimum amplitude threshold")
    parser.add_argument("--t-threshold", type=float, default=0.05, help="t-test p-value threshold")
    parser.add_argument("--binom-threshold", type=float, default=0.0625, help="Binomial test p-value threshold")
    parser.add_argument("--show-plots", action="store_true", help="Display plots interactively")
    parser.add_argument("--no-save-plots", action="store_true", help="Don't save plots to disk")
    args = parser.parse_args()
    
    analyzer = RotationAnalyzer(data_file=args.data, output_dir=args.output)
    results = analyzer.analyze_all_sources(
        sources=args.sources,
        p0=args.p0,
        diff_threshold=args.min_amp,
        t_test_threshold=args.t_threshold,
        binom_threshold=args.binom_threshold,
        save_plots=not args.no_save_plots,
        show_plots=args.show_plots
    )
    
    summary = analyzer.get_summary_statistics()
    print("\n" + "="*50)
    print("ANALYSIS SUMMARY")
    print("="*50)
    for key, value in summary.items():
        if isinstance(value, dict):
            print(f"{key}:")
            for sub_key, sub_value in value.items():
                print(f"  {sub_key}: {sub_value}")
        else:
            print(f"{key}: {value}")

if __name__ == "__main__":
    main()

# """
# Main analysis module for EVPA rotation detection.
# """

# import os
# import numpy as np
# import pandas as pd
# from scipy import stats
# from tqdm import tqdm

# from .loader import load_data, prepare_source_data
# from .bayesian_blocks import (
#     compute_bayesian_blocks, 
#     compute_local_extrema,
#     get_bin_boundaries_for_interval
# )
# from .plotter import plot_rotation_event


# class RotationAnalyzer:
#     """
#     Main class for analyzing EVPA rotations in polarimetric monitoring data.
#     """
    
#     def __init__(self, data_file=None, data=None, output_dir="rotation_analysis"):
#         """
#         Initialize the rotation analyzer.
        
#         Parameters:
#         -----------
#         data_file : str, optional
#             Path to CSV data file
#         data : pd.DataFrame, optional
#             Pre-loaded data
#         output_dir : str, default="rotation_analysis"
#             Directory to save output files and plots
#         """
#         if data is not None:
#             self.data = data
#         elif data_file is not None:
#             self.data = load_data(data_file)
#         else:
#             raise ValueError("Either data_file or data must be provided")
        
#         self.output_dir = output_dir
#         os.makedirs(output_dir, exist_ok=True)
#         os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)
        
#         # Storage for results
#         self.rotation_data = []
#         self.recorded_rotations = set()
        
#     def analyze_source(self, source_name, p0=0.001, diff_threshold=0, 
#                       t_test_threshold=0.05, binom_threshold=0.0625,
#                       save_plots=True, show_plots=False):
#         """
#         Analyze a single source for rotation events.
        
#         Parameters:
#         -----------
#         source_name : str
#             Name of the source to analyze
#         p0 : float, default=0.001
#             Prior probability for Bayesian Blocks
#         diff_threshold : float, default=0
#             Minimum EVPA amplitude to consider as rotation
#         t_test_threshold : float, default=0.05
#             p-value threshold for t-test
#         binom_threshold : float, default=0.0625
#             p-value threshold for binomial test
#         save_plots : bool, default=True
#             Whether to save rotation plots
#         show_plots : bool, default=False
#             Whether to display plots
            
#         Returns:
#         --------
#         list
#             List of detected rotation events for this source
#         """
#         print(f"\nAnalyzing source: {source_name}")
        
#         df_source = prepare_source_data(self.data, source_name)
#         total_periods = df_source['segment'].nunique()
        
#         print(f"Found {total_periods} observation periods")
        
#         source_rotations = []
        
#         for period_index in df_source['segment'].unique():
#             period = df_source[df_source['segment'] == period_index].copy()
            
#             if len(period) < 5:  # Need minimum points for meaningful analysis
#                 continue
                
#             print(f"  Period {period_index}: {len(period)} points")
            
#             period_rotations = self._analyze_period(
#                 period, source_name, period_index, p0, diff_threshold,
#                 t_test_threshold, binom_threshold, save_plots, show_plots
#             )
            
#             source_rotations.extend(period_rotations)
        
#         return source_rotations
    
#     def _analyze_period(self, period, source_name, period_index, p0, diff_threshold,
#                        t_test_threshold, binom_threshold, save_plots, show_plots):
#         """
#         Analyze a single observation period for rotations.
#         """
#         period_rotations = []
        
#         # Remove duplicates
#         period = period.drop_duplicates(subset='MJD').sort_values('MJD')
        
#         # Compute Bayesian blocks for the entire period
#         try:
#             edges, bin_means, bin_errors = compute_bayesian_blocks(period, p0=p0)
#         except Exception as e:
#             print(f"    Bayesian blocks failed: {e}")
#             return period_rotations
        
#         # Convert to numpy arrays
#         bin_means = np.array(bin_means)
#         bin_errors = np.array(bin_errors)
#         bin_centers = 0.5 * (edges[:-1] + edges[1:])
        
#         # Find local extrema
#         maxima_dict, minima_dict = compute_local_extrema(bin_means, bin_centers)
        
#         # Combine extrema and sort by time
#         all_extrema = [(m, v, 'max') for m, v in zip(maxima_dict['MJD'], maxima_dict['EVPA'])]
#         all_extrema += [(m, v, 'min') for m, v in zip(minima_dict['MJD'], minima_dict['EVPA'])]
#         all_extrema.sort(key=lambda x: x[0])
        
#         if len(all_extrema) < 2:
#             print(f"    Not enough extrema found ({len(all_extrema)})")
#             return period_rotations
        
#         # Check for overlaps with existing intervals
#         existing_intervals = [
#             (r['t_start'], r['t_end'])
#             for r in self.rotation_data
#             if r['Source'] == source_name
#         ]
        
#         print(f"    Analyzing {len(all_extrema)-1} adjacent extrema pairs")
        
#         # Analyze adjacent extrema pairs
#         for i in range(len(all_extrema) - 1):
#             t_start, evpa_start, type_start = all_extrema[i]
#             t_end, evpa_end, type_end = all_extrema[i + 1]
            
#             # Only analyze max→min or min→max transitions
#             if not ((type_start == 'max' and type_end == 'min') or
#                     (type_start == 'min' and type_end == 'max')):
#                 continue
            
#             # Check for overlap with existing intervals
#             if any(not (t_end <= s or t_start >= e) for s, e in existing_intervals):
#                 continue
            
#             # Initial EVPA difference check
#             diff_evpa = abs(evpa_end - evpa_start)
#             if diff_evpa < diff_threshold:
#                 continue
            
#             # Get block boundaries
#             start_edge, end_edge = get_bin_boundaries_for_interval(edges, t_start, t_end)
            
#             # Extract sub-data for this interval
#             sub_data = period[
#                 (period['MJD'] >= start_edge) & (period['MJD'] <= end_edge)
#             ].sort_values('MJD')
            
#             if len(sub_data) < 4:
#                 continue
            
#             # Compute sub-period Bayesian blocks
#             try:
#                 edges_sub, means_sub, errs_sub = compute_bayesian_blocks(sub_data, p0=p0)
#             except Exception as e:
#                 print(f"      Sub-BB failed for {t_start:.1f}-{t_end:.1f}: {e}")
#                 continue
            
#             if (len(edges_sub) - 1) <= 2:
#                 continue
            
#             # Statistical tests
#             t_stat, t_pval = stats.ttest_1samp(sub_data['adjusted_EVPA'], popmean=180)
#             n_above = (sub_data['adjusted_EVPA'] > 180).sum()
#             binom_p = stats.binomtest(n_above, len(sub_data), p=0.5).pvalue
            
#             # Amplitude calculation
#             amp = sub_data['adjusted_EVPA'].max() - sub_data['adjusted_EVPA'].min()
#             emax = sub_data.loc[sub_data['adjusted_EVPA'].idxmax(), 'err_EVPA[deg]']
#             emin = sub_data.loc[sub_data['adjusted_EVPA'].idxmin(), 'err_EVPA[deg]']
#             amp_err = np.hypot(emax, emin)
            
#             if amp < diff_threshold:
#                 continue
            
#             # Check if this qualifies as a rotation
#             rotation_key = (source_name, round(start_edge, 1), round(end_edge, 1))
            
#             if (t_pval <= t_test_threshold and 
#                 binom_p <= binom_threshold and 
#                 amp >= diff_threshold and 
#                 rotation_key not in self.recorded_rotations):
                
#                 # Record the rotation
#                 rotation_event = {
#                     'Source': source_name,
#                     'Period': period_index,
#                     't_start': start_edge,
#                     't_end': end_edge,
#                     'ΔMJD': end_edge - start_edge,
#                     'AMPLITUDE': amp,
#                     'AMPLITUDE_error': amp_err,
#                     't-test p': t_pval,
#                     'Binom test p': binom_p,
#                     'n_points': len(sub_data),
#                     'extrema_type': f"{type_start}→{type_end}"
#                 }
                
#                 self.rotation_data.append(rotation_event)
#                 self.recorded_rotations.add(rotation_key)
#                 existing_intervals.append((start_edge, end_edge))
#                 period_rotations.append(rotation_event)
                
#                 print(f"      ROTATION DETECTED: MJD {start_edge:.1f}-{end_edge:.1f}, "
#                       f"Amp={amp:.1f}°, t-p={t_pval:.2e}, binom-p={binom_p:.2e}")
                
#                 # Create plot
#                 if save_plots or show_plots:
#                     plot_output_dir = os.path.join(self.output_dir, "plots") if save_plots else None
                    
#                     plot_rotation_event(
#                         sub_data, period, source_name, period_index,
#                         start_edge, end_edge, amp, t_pval, binom_p,
#                         edges, bin_means, bin_errors, plot_output_dir
#                     )
                    
#                     if not show_plots:
#                         import matplotlib.pyplot as plt
#                         plt.close()
        
#         return period_rotations
    
#     def analyze_all_sources(self, sources=None, **kwargs):
#         """
#         Analyze all sources (or specified subset) for rotation events.
        
#         Parameters:
#         -----------
#         sources : list, optional
#             List of source names to analyze. If None, analyze all sources.
#         **kwargs : dict
#             Additional parameters passed to analyze_source()
            
#         Returns:
#         --------
#         pd.DataFrame
#             Summary of all detected rotation events
#         """
#         if sources is None:
#             sources = self.data['J2000_name'].unique()
        
#         print(f"Starting analysis of {len(sources)} sources...")
        
#         all_rotations = []
        
#         for source in tqdm(sources, desc="Analyzing sources"):
#             try:
#                 source_rotations = self.analyze_source(source, **kwargs)
#                 all_rotations.extend(source_rotations)
#             except Exception as e:
#                 print(f"Error analyzing {source}: {e}")
#                 continue
        
#         print(f"\nAnalysis complete. Found {len(all_rotations)} rotation events.")
        
#         # Save results
#         if self.rotation_data:
#             results_df = pd.DataFrame(self.rotation_data)
#             output_file = os.path.join(self.output_dir, "rotation_events.csv")
#             results_df.to_csv(output_file, index=False)
#             print(f"Results saved to: {output_file}")
#             return results_df
#         else:
#             print("No rotation events detected.")
#             return pd.DataFrame()
    
#     def get_summary_statistics(self):
#         """
#         Generate summary statistics of detected rotations.
        
#         Returns:
#         --------
#         dict
#             Dictionary with summary statistics
#         """
#         if not self.rotation_data:
#             return {"message": "No rotation events detected"}
        
#         df = pd.DataFrame(self.rotation_data)
        
#         summary = {
#             "total_rotations": len(df),
#             "unique_sources": df['Source'].nunique(),
#             "sources_with_rotations": df['Source'].unique().tolist(),
#             "amplitude_stats": {
#                 "mean": df['AMPLITUDE'].mean(),
#                 "std": df['AMPLITUDE'].std(),
#                 "min": df['AMPLITUDE'].min(),
#                 "max": df['AMPLITUDE'].max(),
#                 "median": df['AMPLITUDE'].median()
#             },
#             "duration_stats": {
#                 "mean": df['ΔMJD'].mean(),
#                 "std": df['ΔMJD'].std(),
#                 "min": df['ΔMJD'].min(),
#                 "max": df['ΔMJD'].max(),
#                 "median": df['ΔMJD'].median()
#             },
#             "rotations_per_source": df['Source'].value_counts().to_dict()
#         }
        
#         return summary


# def main():
#     """
#     Main function to run the analysis from command line.
#     """
#     import argparse
    
#     parser = argparse.ArgumentParser(description="EVPA Rotation Analysis")
#     parser.add_argument("--data", required=True, help="Path to CSV data file")
#     parser.add_argument("--output", default="rotation_analysis", 
#                        help="Output directory")
#     parser.add_argument("--sources", nargs="+", 
#                        help="Specific sources to analyze (default: all)")
#     parser.add_argument("--p0", type=float, default=0.001,
#                        help="Prior probability for Bayesian blocks")
#     parser.add_argument("--min-amp", type=float, default=0,
#                        help="Minimum amplitude threshold")
#     parser.add_argument("--t-threshold", type=float, default=0.05,
#                        help="t-test p-value threshold")
#     parser.add_argument("--binom-threshold", type=float, default=0.0625,
#                        help="Binomial test p-value threshold")
#     parser.add_argument("--show-plots", action="store_true",
#                        help="Display plots interactively")
#     parser.add_argument("--no-save-plots", action="store_true",
#                        help="Don't save plots to disk")
    
#     args = parser.parse_args()
    
#     # Initialize analyzer
#     analyzer = RotationAnalyzer(data_file=args.data, output_dir=args.output)
    
#     # Run analysis
#     results = analyzer.analyze_all_sources(
#         sources=args.sources,
#         p0=args.p0,
#         diff_threshold=args.min_amp,
#         t_test_threshold=args.t_threshold,
#         binom_threshold=args.binom_threshold,
#         save_plots=not args.no_save_plots,
#         show_plots=args.show_plots
#     )
    
#     # Print summary
#     summary = analyzer.get_summary_statistics()
#     print("\n" + "="*50)
#     print("ANALYSIS SUMMARY")
#     print("="*50)
#     for key, value in summary.items():
#         if isinstance(value, dict):
#             print(f"{key}:")
#             for sub_key, sub_value in value.items():
#                 print(f"  {sub_key}: {sub_value}")
#         else:
#             print(f"{key}: {value}")


# if __name__ == "__main__":
#     main()

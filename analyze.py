from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(data_file="data/monitoring_data.csv")
results = analyzer.analyze_all_sources(
    p0=0.001, diff_threshold=90,
    t_test_threshold=0.05, binom_threshold=0.0625
)
print(analyzer.get_summary_statistics())

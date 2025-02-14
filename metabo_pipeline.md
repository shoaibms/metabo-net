📦 Metabolomics-Analysis-Pipeline
 ├── 📂 data                          
 │
 ├── 📂 src                           
 │   ├── 📂 1_data_preprocessing    # Data pre-processing scripts
 │   │   ├── feature_filter.py       # Initial feature filtering
 │   │   ├── missing_viz.py          # Missing value visualization
 │   │   ├── mar_test.py             # Missing at Random test
 │   │   ├── logistic_test.py        # Logistic regression test
 │   │   ├── mcar_test.py            # Missing Completely at Random test
 │   │   ├── median_impute.py        # Median imputation
 │   │   ├── rf_impute.R             # Random Forest imputation
 │   │   ├── ml_impute.py            # Machine learning imputation
 │   │   ├── impute_validate.py      # Imputation validation
 │   │   ├── impute_dist_check.py    # Distribution check after imputation
 │   │   ├── isolation_forest.py      # Isolation Forest for outliers
 │   │   ├── dim_reduce_outliers.py   # Dimensionality reduction for outliers
 │   │   ├── outlier_viz.py           # Outlier visualization
 │   │   ├── transform_data.py        # Data transformation
 │   │   ├── normality_test.py        # Normality testing
 │   │   ├── normality_viz.py         # Normality visualization
 │   │   ├── transform_metrics.py      # Transformation metrics
 │   │   ├── transform_eval.py         # Transformation evaluation
 │   │   ├── variance_calc.py          # Variance calculation
 │   │   └── diversity_metrics.py      # Diversity metrics calculation
 │   │
 │   ├── 📂 2_analysis                # Main analysis scripts
 │   │   ├── ini_analysis.py            # Final preprocessing
 │   │   ├── ini_analysis_summary.py    # Preprocessing summary
 │   │   ├── stat_tests.py            # Statistical tests
 │   │   ├── pls_tissue.py            # PLS analysis by tissue
 │   │   ├── spearman_network.py      # Spearman correlation network
 │   │   ├── network_decay.py         # Network decay analysis
 │   │   ├── network_summary.py       # Network summary
 │   │   ├── baysian_network.R        # baysian network
 │   │   ├── tissue_analysis.R        # Tissue-specific analysis
 │   │   └── tissue_summary.R         # Tissue analysis summary
 │   │
 │   ├── 📂 3_visualization           # Plotting scripts
 │   │   ├── 📂 figure1              
 │   │   │   ├── network_viz.py      # Network visualization
 │   │   │   ├── radar_plot.R        # Radar plot
 │   │   │   ├── bayesian_crosstalk.R # Bayesian network analysis
 │   │   │   ├── hub_dist.R          # Hub distribution
 │   │   │   ├── hub_decay.R         # Hub decay
 │   │   │   ├── module_org.R        # Module organization
 │   │   │   ├── module_stability.R  # Module stability
 │   │   │   └── temporal_stability.R # Temporal stability
 │   │   │
 │   │   ├── 📂 figure2
 │   │   │   ├── tissue_temporal.R   # Tissue temporal analysis
 │   │   │   ├── temporal_corr.R     # Temporal correlation
 │   │   │   └── tissue_plot.R       # Tissue plotting
 │   │   │
 │   │   └── 📂 figure3
 │   │       └── validation_viz.R     # Validation visualization
 │   │
 │   └── 📂 4_chemical_identification # Chemical ID scripts
 │       ├── hmdb_annotate.py         # HMDB annotation
 │       ├── gnps_annotate.py         # GNPS annotation
 │       ├── struct_classify.py       # Structural classification
 │       └── func_group.py            # Functional group analysis
 │
 ├── 📂 3D_figures                    # Interactive plot
 ├── 📂 docs                          # Documentation 
 ├── requirements.txt                  # Dependencies
 ├── environment.yaml                  # Configuration 
 └── README.md                         # Project overview


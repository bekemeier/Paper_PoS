The paper "Cerebral Encoding of Parts of Speech is Distributed and Context-Dependent" employs the publicly available MOUS dataset (Schoffelen et al., 2019) (https://doi.org/10.34973/37n0-yc51) and the analysis pipeline implemented by Huizeling et al. (2021) (https://doi.org/10.1162/nol_a_00054). In contrast to the latter study, we include the factor Part-of-Speech (POS) into our analyses. 
The preprocessing, LCMV beamforming, and multiset canonical correlation analysis (MCCA) follow the analysis pipeline published by Arana et al. (2020; https://doi.org/10.34973/tf8r-rq72). Specifically, we follow their pipeline up to the MCCA stage (implemented in mous_supramodal_JoN_pipeline.m), after which we perform encoding-model analyses. The subsequent analysis steps (Steps 2–10 below) are implemented in our own scripts and conceptually follow the encoding-model framework described by Huizeling et al. (2021).
The analysis pipeline is as follows:
1.	"mous_supramodal_JoN_pipeline.m" (original pipeline by Arana et al. (2020): (https://doi.org/10.34973/tf8r-rq72)): 
		(i) retrieve which subject belongs to which scenario, load subject data for the respective modality; 
		(ii) Compute parcelwise source time-course (LCMV + PCA); 
		(iii) Perform multiset CCA to obtain canonical components and save the resulting tlck structure, containing source-level activity time-locked to individual word onsets together with the corresponding trial information.

2. 	"CrossvalRidge_MOUS_WCl.m":
		(i) the script fits cross-validated Ridge Regression models for the original data and for the shuffled data (50 permutations - to generate a distribution under the null hypothesis).
		(ii) saves prediction error and regression coefficients for (1) the original data, (2) the shuffled data, and (3) averaged error and coefficients for the shuffled data.
		(iii) This procedure is performed separately for 
		- the full model (all predictors), and
		- the reduced model (all predictors excluding the predictor of interest).
		
3.	"ModelComparison_MOUS.m":
		This script computes the coefficient of determination (COD or R^2), defined as: COD = 1-(errorFullModel/errorReducedModel.
		It first loads the error matrices (time x subject) for both full and reduced models and computes individual COD values per timepoint and subject, which are saved for statistical analysis. 
		In addition, COD values are computed for the corresponding shuffled (null) models using the same procedure, producing a null distribution for statistical inference.

4.	"combineSce_forstats.m":
		The script aggregates COD values across all scenarios, combining both original and shuffled data into a single dataset. 
		The resulting matrices have dimensions time X subject.
5.	"MOUS_stats.m":
		Nonparametric permutation tests using a dependent-samples t-statistic across subjects. 
		Multiple comparisons across time are accounted for by using a max-statistic procedure based on 5000 permutations.
6.	"FDR_post_perm.m":
		Spatial multiple-comparisons correction across cortical parcels using the Benjamini–Hochberg FDR procedure.
		The script applies FDR correction to p-values derived from the time-corrected permutation statistics, controlling for multiple comparisons across space (parcels).
		Significant parcels are identified as those containing at least one timepoint surviving both temporal permutation correction and spatial FDR correction.
7.	"Plot_MOUS_loop.m":
		This script visualizes t-statistics on cortical surfaces separately for each hemisphere and for lateral and medial views.
		Maps are thresholded using the minimum t-value that survives temporal permutation testing and spatial FDR correction (t_thresh). Significant values are projected onto cortical surface maps for each hemisphere using a colormap.

8. 	"post_hocs_MOUS.m"
		The script performs time-resolved non-parametric paired permutation tests on regression coefficients extracted from significant parcels. 
		For each timepoint, a paired t-statistic is computed across subjects comparing two conditions (e.g., verb vs adjective, verb vs noun, adjective vs noun). 
		Multiple comparisons across time are accounted for by using a max-statistic procedure based on 5000 permutations.
 		Significant timepoints are identified at α = 0.0167 (Bonferroni-corrected for three comparisons).
9. 	"FDR_correction_posthocs.m"
		This script aggregates time-resolved post-hoc permutation p-values across all previously identified significant parcels and performs spatial multiple-comparisons correction. 
		For each parcel, the minimum uncorrected p-value across time is extracted for each contrast (verb vs noun, adjective vs noun, verb vs adjective), yielding one summary statistic per parcel. 
		These parcel-level p-values are then corrected for multiple comparisons across space using the Benjamini–Hochberg false discovery rate (FDR) procedure. 
		Parcels surviving FDR correction are identified as showing reliable post-hoc effects, and their corresponding significant timepoints (uncorrected p < 0.0167) are retained for interpretation and visualization.
10. 	"Plot_PostHoc_curves.m"
		This script reconstructs parcel-wise post-hoc PoS effects and visualizes their temporal evolution as waveforms, allowing inspection of when and how specific PoS contrasts emerge over time in individual cortical regions.
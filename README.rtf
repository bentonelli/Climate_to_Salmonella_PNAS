{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww19860\viewh14820\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 README for Tonelli et al. 2025\
\
Please note that some data files used to conduct this analysis need to be download from external sources due to copyright or other restrictions. As a result, some intermediate data files created in the analysis retain some original data. We indicate these data files with the following notation:\
\
* Indicates when data needs to be downloaded from external sources. Links provided.\
** Indicates when data needs to be derived by the user from external sources using code provided here\
\
Code:\
\
1_data_processing: Code to download and process raw data into a format for input into the model\
\
	- 1a_download_daymet.R: Download daymet using associated package \'93daymetr\'94\
	- 1b_process_daymet.R: Take daymet data and calculate difference in temperatures for each consecutive July.\
	- 1c_deltaT_to_cell.R: Convert delta-t rasters into cell format\
	- 1d_er_by_cell.R: Classify cells as belonging to a given ecoregion\
	- 1e_process_tree_cover.R: Using a pre-constructed mosaic raster of tree cover at the continental level, calculate the % tree cover of each spatial cell\
	- 1f_process_irruptions.R: Converts Christmas Bird Count data into irruption time series for 8 species across 2 regions\
	- 1g_irruption_historical.R: Does the same as 1f, but creates time series blinded to certain years \
	- 1h_process_whispers.R: Convert WHISPers data into regional time series of the number of individual animals reported in salmonellosis outbreaks \
	- 1i_compile_mdl_data.R: Formats cone production data, then takes time series information, and packages it to be read into a Stan model.\
	- 1j_trim_model_data.R: Does the same as 1i, but creates historical data for the forecast evaluation.\
\
2_model_run: Code to run hypothesis and historical models\
\
	- 2a_cmid.stan: Bayesian model written in Stan, used for both hypothesis-testing and historical models\
	- 2b_run_model.R: Code to run hypothesis-testing model\
	- 2c_run_time_maching_model.R: Code to run historical models\
 \
3_model_analysis: Code to analyze the model output, and parts of the data used for the model\
\
	- 3a_mdl_main_analysis.R: Extracts parameter estimates for those reported in the manuscript\
	- 3b_effect_size_analysis.R: Extract effect sizes for certain parameters of interest.\
	- 3c_PPCs.R: See how well model reproduces cone production data.\
	- 3d_mdl_check.R: See if model passes diagnostic checks.\
	- 3e_PPO_check.R: Run posterior predictive overlap checks for hypothesis-testing model\
	 \
	- Figures: Folder containing scripts to recreate figures in the manuscript\
		- 3e_Fig2_main_results.R\
		- 3f_Fig3_disease_data.R\
		- 3g_FigS1_methodological_overview.R\
		- 3h_FigS2_cat_plots.R\
		- 3i_FigS5_data_availability.R\
\
4_predictive_analysis: Code to score forecasts based on accuracy, and evaluate these predictions visually\
\
	-4a_predictive_scoring.R: Reads in all historical models to score probabilistic forecasts.\
	-4b_forecasting_plot.R: Calculate irruption likelihood scores for public-facing forecasts.\
	-4c_brier_score.R: Calculate brier scores.\
\
Data:\
\
1_raw: Contains raw data downloadable from external sources\
	- CBC - CA and US data provided by request from Audubon Society. https://netapp.audubon.org/cbcobservation/\
		- *US_CA_CBC_data.csv\
	- Counties - Data from https://www.weather.gov/gis/Counties\
		- *c_05mr24.shp\
	- Daymet - Downloaded 10/24 via R package \'93daymetr\'94 with \'93download_daymet.R located in code/1_data_processing\
		- *Files in format: tmax_july_[Year]_ncss.nc\
	- Ecoregions - US EPA ecoregion shape file, dowloaded via: https://www.epa.gov/eco-research/ecoregions-north-america\
		- *NA_CEC_Eco_Level1.shp\
	- MASTREE - Downloaded 10/24 from Github repo: https://github.com/JJFoest/MASTREEplus/blob/main/Data/MASTREEplus_2024-06-26_V2.csv\
		- *MASTREEplus_2024-06-26_V2\
	- Other - Contains one file that encodes whether each animal mentioned in the WHISpers database is identified to species, is a passerine, and is a bird.\
		- animal_species.csv\
	- Tree_cover - Contains one file in .nc format that encodes the % forest cover at 1km resolution from the year 2000. Downloaded in tile format from: https://lpdaac.usgs.gov/products/gfcc30tcv003. After download at 30m resolution, data tiles were merged and aggregated to 1km resolution. Note that there is one small gap in data availability in central Alaska. \
		- *FC_2000_1km.nc\
	- WHISPers - Downloaded 10/24 from https://whispers.usgs.gov/home. All events with Salmonellosis as a possible diagnosis between Nov.1, 1987 and Oct. 1 2024.\
		- *whispers_raw.csv\
2_formatted: Contains formatted data for DeltaT, Disease, Ecoregion, and Tree Cover.\
	- DeltaT\
		- deltaT_cells: Files in .rds format encoding the temperature difference between each pair of consecutive years, averaged by cell identity.\
		- **deltaT_rasters: Files in .nc format encoding the temperature difference between each pair of consecutive years, as rasters.\
	- Disease\
		- **disease_locations.rds: Disease event information from WHISPers with location information retained.\
		- east_salmonellosis_ts.rds: Yearly time series of the number of individuals involved in salmonellosis outbreaks in the eastern region. \
		- spec_freq_e.csv: The number of times each species was involved in suspected salmonellosis outbreaks, eastern region.\
		- spec_freq_w.csv: The number of times each species was involved in suspected salmonellosis outbreaks, eastern region.\
		- west_salmonellosis_ts.rds: Yearly time series of the number of individuals involved in salmonellosis outbreaks in the western region. \
	- Ecoregions\
		- ER_by_cell.rds: Cell ecoregion identity at the level for masting level of the analysis.\
		- ER_by_cell8.rds: Cell ecoregion identity at the level for irruption level of the analysis.\
	- Irruptions\
		- **cbc_core_w_cell_filt.rds: Processed, intermediate CBC data used to construct regional time series of irruptions.\
		- irruptions_detrended_east.csv: Processed irruption time series for eastern species.\
		- irruptions_detrended_west.csv: Processed irruption time series for western species.\
		- Time_machine: Processed irruption time series used for forecasting model testing.\
	- Tree_cover\
		- FC_by_cell.rds: Forest cover by cell identity.\
\
3_model_in: Data to be read into Bayesian statistical models\
	- model_full_data_2024.rds: Model data used for the hypothesis testing model.\
	- time_machine: Model data formatted in the same way as the hypothesis testing model, but based on historical datasets. Available on request due to large file size (>50GB)\
		- blinded: Model data blinded to contemporary tree masting information.\
		- no_blinded: Model data not blinded to contemporary tree masting information.\
\
4_model_out: Fit model objects from hypothesis testing model, and historical runs.\
	- mdl_fit.rds: Fit model object from hypothesis testing model.\
	- time_machine: Contains fit model objects for blinded and unblinded historical models.\
\
5_predictive_analysis\
	- dis_prediction_scoring_blinded.rds - Scoring for blinded disease predictions\
	- dis_prediction_scoring_cascade.rds - Scoring for non blinded disease predictions\
	- irr_prediction_scoring_blinded.rds - Scoring for blinded irruptions predictions\
	- irr_prediction_scoring_cascade.rds - Scoring for non blinded irruptions predictions\
	- pred_analysis_blinded.pdf - Individual year prediction graphics, blinded time machine models\
 	- pred_analysis_not_blinded.pdf - Individual year prediction graphics, non-blinded time machine models\
\
\
\
}
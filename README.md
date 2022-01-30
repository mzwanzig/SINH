# SINH
An extended Michaelis-Menten model that includes the effect of substrate inhibition (SINH)
* presenting model definition
* exemplary usage
* data demonstrating subtrate inhibition effects for 3 different enzymes (-> kinetics of 3 x 120 deadwood samples)
* reproduce findings published in this manuscript:

Zwanzig, M., Nähte, K., Michalzik, B., Tischer, A., An extended Michaelis-Menten model that includes substrate inhibition for improved kinetic description of extracellular hydrolytic enzymes in deadwood. (Submitted to Soil Biology and Biochemistry 2022)
* please check this manuscript for more information on the model behavior and the interpretation of its results

## /data
* contains a single file "pk1_hydrolases_sep.xlsx" representing the enzyme-substrate reaction data of 120 samples for each of the three EHE presented in the article mentioned above

## /output
* contains the results of the analysis done by "fit_and_evaluate_nls_EHE_models.R" (see below)

## /images
* contains R-scripts creating the subfigures shown in the article mentioned above
* in addition, previews of these figures are given as pdf files
* (R-scripts depend on files in /output)

## fit_and_evaluate_nls_EHE_models.R
i) loads the complete data
ii) creates a subset representing the substrate enzyme kinetics of one type of enzyme for a single sample (including 3 replicates per substrate concentration tested)
iii) fits the classical Michaelis Menten model (MM)
iv) fits the substrate inhibition model (INH)
v) extracts parameter estimates and evaluates the performance of each type of model
vi) creates a figure showing each model fit and its respective residuals in /output
vii) step ii-vi are repeated for all other samples and enzymes
viii) writes all results of step iv) in separate csv-files in /output

## example.R
* uses a small exemplary dataset to demonstrate the fitting of the substrate inhibition model SINH and its comparison to the classical Michaelis-Menten approach
* to analyze other data, replace the data object by one of your relating a vector of substrate concentrations (named "S") to a vector containing reaction rates (named "R")

## example.pdf
* comparison of MM and SINH as given by "example.R" and "plot_MMandSINH.R"

## plot_MMandSINH.R
* routine to fit, bootstrap, calculate and plot both the classical Michaelis-Menten (MM) and the substrate inhibition model (SINH) for a given data.frame named "df_i" that contains a column "S" (substrate concentration) and "R" (reaction rate)
* used by "example.R" and "images/Figure1.R"

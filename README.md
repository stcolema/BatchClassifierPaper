# BatchClassifierPaper
Repository for the analysis of [Coleman et al., 2022](https://doi.org/10.1101/2022.01.14.476352). Note that all scripts assume that the ``BatchClassifierPaper`` directory is the working directory.

Many of the scripts in this repo require some of the lead author's R package. To install these please run the following lines of code:

```
devtools::install_github("stcolema/BatchMixtureModel", ref = "Paper")
devtools::install_github("stcolema/mdiHelpR")
```

## Simulation study
To recreate the simulation study analysis, the scripts that are in the ``Scripts/Simulations/`` directory are called in the following order:

1. ``dataGeneration.R``: generate ten datasets in each of the six scenarios and save the scenario description file.
2. ``modelling.R``: run the models and save various outputs.
3. ``acceptance.R``: check that acceptance rates are close to the [0.1, 0.5] range (good target for efficient exploration).
4. ``gewekeConvergence.R``: use the Geweke statistic for the complete log-likelihood to check within-chain convergence.
5. ``completeLikelihood.R``: use the complete log-likelihood to check across-chain convergence. This is a manual step and we identify chains that have achieved within-chain convergence but have settled in a different mode to the other chains in the same simulation and any chains that incorrectly passed the previous step.
6. ``modelComparison.R``: visualise the performance of the various models having dropped the poorly-behaved chains.
7. ``lookAtGeneratedData.R``: this produces a plot of the observed and batch-corrected datasets from one example for each scenario. These are used in the Supplementary Material.

## Dopico et al. analysis
To recreate the analysis of the ELISA data from Dopico et al. (2021), please download the dataset ``https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData`` to ``./Data/ELISA/Sweden/``. Then run 

1. ``swedenModelling.R``: run MCMC for the MVT mixture model.
2. ``swedenModelCheck.R``: assess chains on log-likelihood trace plots.
3. ``swedenHyperParameterAndSeroprevalencePlot.R``: visualise the seroprevalence estimate and the prior distributions.
4. ``swedenDataPlots.R``: compare the observed data to the inferred.

## ELISA like simulations
To recreate the simulation study of data generated from the MVT model applied to the ELISA data from Dopico et al. (2021),  the scripts that are in the ``Scripts/pseudo-ELISA/`` directory are called in the following order

1. ``elisaSimGen.R``: generate the data.
2. ``elisaLikeModelling.R``: apply each of the models and save various outputs.
3. ``elisaLikeConvergence.R``: check within-chain convergence using the Geweke-statistic.
4. ``elisaLikeLikelihood.R``: plot the sampled likelihoods and manually remove chains that found a local minimim.
5. ``elisaLikeModelComparison.R``: compare the model performance under the F1 score and the distance.
6. ``elisaFinalResults.R``: make the plots actually used in the paper.

## Dingens et al. analysis
To recreate the analysis of the ELISA data from Dingens et al. (2020), please run 

1. ``seattleModelling.R``: run multiple MCMC chains on the data.
2. ``seattleModelCheck.R``: assess convergence using the log-likelihood.
3. ``seattleResultsPlot.R``: create the figure used in the paper.

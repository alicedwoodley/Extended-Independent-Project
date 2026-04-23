<p align="center">
𓆝 𓆟 𓆞 𓆟 𓆝
</p>

# MMath Extended Independent Project (in progress)

### Overview
In this project, I propose multispecies stock-recruitment models that directly incorporate the effects of predator-prey interactions. I am producing a written report and will additionally give a talk on the results of my work.

### Data
This project uses data from the RAM Legacy stock assessment database (<https://www.ramlegacy.org/>). The folder 'RAM Legacy raw data' contains both the Rdata file and script file describing the contents of the database. 

## Repository contents

### 0 - Initial research
Here I leave some of my workings following chapter 13 of '[Introductory Fisheries Analyses with R](https://derekogle.com/IFAR/)' by Derek Ogle. I made basic plots and fitted single-species stock-recruitment models, first for the Pacific Coast arrowtooth flounder and then for the North Sea European plaice. \
(Note: this code is not discussed in my final report, but forms the basis for fitting the included stock-recruitment models.)

### 1 - Selecting species
Here I search the database for a suitable group of species whose data I can use to develop models, settling on the Gulf of Alaska. 

### 1a - Data processing
The data of the chosen species from the Gulf of Alaska is processed here, in preparation for use in model development. All further work uses the RData file _Gulf_of_Alaska_ created at the end of processing and included in this repository.

### 2a - Single species Ricker models
Traditional Ricker models are fit here to investigate how prominent a stock-recruitment relationship is in each species.

### 2b - Multispecies Ricker models
Here I experiment with adapting a Ricker model to incorporate SSB of interacting species as additional explanatory variables and compare these new models to the single species alternatives. \
(Note: not all of these models are included in the final report.)

### Additional - Multilevel linear regression
Here a multilevel linear regression model was used to incorporate all 6 species in a single model. All 6 regression lines had very similar gradient- it is unclear as to why that is. This model could be improved to incorporate a covariate that represents how much each species is preyed on / how much of their food source is available. \
(Note: this model was one of my first experimental ideas and does not feature in the final report.)

### Report and presentation
The final report will be complete by 27/04/2026 and this work will be presented on 08/05/2026. Both the completed report and presentation slides will be uploaded here once complete.

<p align="center">
𓆝 𓆟 𓆞 𓆟 𓆝
</p>

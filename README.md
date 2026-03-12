<p align="center">
𓆝 𓆟 𓆞 𓆟 𓆝
</p>

# MMath Extended Independent Project (in progress)

### Overview
In this project I propose multispecies stock-recruitment models that directly incorporate the effects of interspecies interactions. I will then produce a written report and additionally give a talk on the results of my work.

### Data
This project uses data from the RAM Legacy stock assessment database (<https://www.ramlegacy.org/>). The folder 'RAM Legacy raw data' contains both the Rdata file and script file describing the contents of the database. 

## Repository contents

### 0 - Initial research
Here I leave some of my workings following chapter 13 of '[Introductory Fisheries Analyses with R](https://derekogle.com/IFAR/)' by Derek Ogle. I made basic plots and fitted single-species stock-recruitment models, first for the Pacific Coast arrowtooth flounder and then for the North Sea European plaice. (Note: this code is not discussed in my final report, but forms the basis for fitting the included stock-recruitment models.)

### 1 - Selecting species
Here I search the database for a suitable group of species whose data I can use to develop models, settling on the Gulf of Alaska. 

### 1a - Data processing
The data of the chosen species from the Gulf of Alaska is processed here, in preparation for use in model development.

### 2a - Multilevel linear regression
Here a multilevel linear regression model was used to incorporate all 6 species in a single model. All 6 regression lines had very similar gradient- it is unclear as to why that is. This model could be improved to incorporate a covariate that represents how much each species is preyed on / how much of their food source is available. (Note: this model may not feature in the final report and if so, may be removed from this repository.)

### 2b - Single species Ricker models
Traditional Ricker models are fit here to investigate how prominent a stock-recruitment relationship is in each species.

### 2c - Multispecies Ricker models
Here I experiment with adapting a Ricker model to incorporate additional explanatory variables, namely the spawning stock biomass of species that interact with walleye pollock and arrowtooth flounder, and compare these new models to the single species alternatives.

### Report and presentation
Report deadline 27/04/2026 and presentation day 08/05/2026. Both the completed report and presentation slides will be uploaded here once complete.

<p align="center">
𓆝 𓆟 𓆞 𓆟 𓆝
</p>

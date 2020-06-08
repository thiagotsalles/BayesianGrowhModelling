# Bayesian growth modelling
Research project aimed at modelling the growht (weight) of animals via frequentist and Bayesian approaches.\
This study is published in Livestock Science.\
All details about methodology and results are on
[https://doi.org/10.1016/j.livsci.2020.104115](https://doi.org/10.1016/j.livsci.2020.104115).\
BayesianGrowhModelling.R contains the execution of the project.

## Highlights
• A way of combining sample and prior information in growth modelling is presented.\
• Growth modelling benefits from Bayesian approach.\
• Prior information improves estimates but its predominance may produce inconsistencies.

## Abstract
Growth models are used to understand the relationships in production during the life of an animal, being an abstraction of their natural dynamics. In this context, the objective of this research was to fit a curve for weight of Santa Ines sheep using frequentist and Bayesian approaches, present strategies for eliciting prior distributions for the latter and compare the results obtained with each one. Growth data from a literature study was used as sample. The parameter estimates were obtained using nonlinear least squares in the frequentist approach and using Monte Carlo method via Markov Chains algorithms in the Bayesian approach. Noninformative and informative prior distributions were used in the Bayesian approach, with prior information coming from other six studies. A methodology for eliciting informative prior distributions was provided. Prior information contributed to more precise estimates of sheep weight. It was seen that predominance of prior information may produce inconsistent interval estimates. Although the values of the parameters estimated by the two approaches were similar, the use of the Bayesian approach, together with the prior distributions, allowed for good and more precise estimates when compared to the frequentist approach.

## Data
The sample comprises recently observed data. The prior information are previous observations of data.\
Sample is used in frequentist approach, while Bayesian approach uses sample combined with prior information.

<img src="./Images/fig1.png" width="600">

## Model
The Gompertz model was used to describe growht:

<img src="./Images/gompertz.png" width="300">

where *y* represents the weight, in kg; *t* is the age, in days; *α*, *β* and *γ* are the model parameters; and *ε* is the random error, being that *ε* ~ *N*(0,σ<sup>2</sup>) .

In Bayesian approach, the prior information was used to define prior distributions of the parameters of the model.
Three options of prior distributions were tested:\
• Non-informative (no influence on model fit).\
• Informative (considerable influence on model fit).\
• Informative with greater variance (moderate influence on model fit).

## Graphical abstract
<img src="./Images/abstract.jpg" width="800">

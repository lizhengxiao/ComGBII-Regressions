# ComGBII-Regressions
This project provide the R codes for the paper '[A new class of composite GBII regression models with varying threshold for modelling heavy-tailed data (arxiv.org)](https://arxiv.org/abs/2203.11469).


- `functions-comGB2-distribution.r` provides the density, distribution function, quantile function and random generation for the composite GBII distribution, along with several risk measures. The initialization methods for the parameters are also included in this file.

- `1-comGB2-model-solnp.r`provides the parameter estimation methods for the composite GBII regression models.

- other `.r` files contains the sub-models dicussed in this paper, including GBIIG, BIIG, BG, IBG, PG and IPG.

- ``Case-I-simulation.r`` shows the simulation study to demonstrate the validity of parameter estimates.

- `comGB2RegM.r` contains estimation results for all seven sub-models proposed in this paper.

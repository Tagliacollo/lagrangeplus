This is an extended version of Richard Ree's Lagrange model (Ree et al 2005). The model estimates rates of species dispersal between areas. Given a phylogeny and the contemporary distribution of the extant species, the model estimates the rates of dispersal between areas. The model also estimates how dispersal rates between areas has changed over time. Estimates are derived using a maximum likelihood approach where the objective is to find the pattern of dispersal that makes the phylogeny and the contemporary distribution of species most likely.

The original Lagrange model used a user-defined matrix of relative dispersal rates and then estimates a scaler value, D, that is multiplied by the relative dispersal matrix to estimate actual dispersal rates. For example, relative dispersal rates can be computed based on the distance between areas and D estimated to give estimates of the true dispersal rates.

We have modified Lagrange so that the terms of a dispersal matrix are estimated independent of each other. This modification allows the model to tell us which connections are most consistent with a given phylogeny, instead of making an assumption of the relative importance of dispersal rates a priori.

We would like to thank Dr. Ree and his colleagues for making the original Lagrange code available. Their openness made the development of our modified model a more tractable task. In the same spirit of openness, we are making our modified version available to the public.

Results based on this model are reported in a manuscript that we are in the process of submitting for publication now. If and when the research is published, we will report the citation here.

References:
Ree, R. H., B. R. Moore, et al. (2005). "A likelihood framework for inferring the evolution of geographic range on phylogenetic trees." Evolution 59(11): 2299-2311.
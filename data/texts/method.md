This app follows the framework proposed in [Herrera Malatesta and de Valeriola (2024)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0307743).
This framework, called the *Robustness Assessing Framework*, focuses on assessing the robustness and uncertainties of conclusions drawn from applying point pattern analysis to, mostly, 
non-systematic regional data in archaeology. 
To achieve this, we have articulated the discussion on the reconstruction of past landscapes using computational methods around three key aspects: 
1) the use of point pattern analysis (PPA) in archaeology
2) the quantification of uncertainties
3) the consideration of robustness.

The framework and subsequently this app are designed to aid archaeologists working with datasets that are known to contain sources of uncertainty in applying
spatial statistical methods and achieve a higher understanding of the uncertainties of the resulting models. 

#### Methodology
The framework itself consists of three simple steps: 

**Step 1 - The Observable**

The first step is defined as "the observable” which is a point clustering metric whose changes are tracked and measured when considering deviations. 
In this case, this is the model resulting from a Pair Correlation Function (PCF) with a Monte Carlo simulation envelope based on 100% of the dataset used as input. 
The app will provide the result of "the observable "which will be used as a reference value for the second step.

**Step 2 - The Experiment**

The second step is "the experiment" where regular intervals of data will be sampled from the loaded dataset. 
More precisely, the app will deduct 10%, 20%, 30%, 40%, and 50% of the database’s sites and perform the PCF again with a Monte 
Carlo simulation envelope for each of the sampled groups. The resulting models are what we have called the "robustness scenarios". 

Note that in our paper we used two sampling methods, one using a uniform distribution and another one using an inhomogeneous distribution. 
To make the app work optimally, only the sampling via the uniform distribution is provided. 
If a more experienced researcher would like to access the code and consider the inhomogeneous distribution, 
please check out the original code <a href="https://osf.io/u2gyq/" target="_blank">here</a>.

**Step 3 - Comparison Tools**

The third step consist of the "comparison tools", which are methods to assess the frequencies and interval midpoint densities. 
The comparison tools is the step that will allow the analyst to assess the robustness and quantify the uncertainty of the spatial models created based on their dataset. 

**Results**

The figure of the first comparison tool will present the percentages of sites that are kept in each robustness scenarios against the percentage of robustness
scenarios in which the conclusion is similar to "the observable". 
This figure provides a direct percentage of the probability that, by extracting a particular percentage of the dataset, 
the results of the robustness scenario will be similar to "the observable".

The second comparison tool goes deeper into understanding these patterns and provides insight into what can be further observed from the clustering patterns.
At this point, it is important to clarify that the app is only analyzing the statistically significant patterns on the PCF, and not the regular ones. 
Again, if a researcher would like to include this in their analysis, they need to modify the source code.

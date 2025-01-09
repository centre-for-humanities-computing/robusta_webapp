On this web app you have two options: upload your own data for processing or use the Montecristo data used by Herrera Malatesta and de Valeriola (2024). The Montecristo data is meant as a test dataset, where you can play a bit around with the parameters and see potential outputs. 

If you would like to work on your own data, then you have to upload the following shapefiles:
1. A vector file of the points to be analysed (e.g., sites, materials). You need to upload a .shp and .shx file.
2. A spatial polygon file of the research area. This also needs to consist of a .shp and .shx file.

Furthermore, the number of Monte Carlo simulations, the number of Robustness Scenarios, and the Quantiles need to be specified. For the Monte Carlo simulations, we recommended values between 1000 and 5000. For the Robustness Scenarios, the value should be between 100 and 3000. Please note that these values are a recommendation based on the standard computational power of a laptop or PC and the possibilities of this web app. But try to play around with it and see what happens. However, if you would like to increase the values to even larger numbers than the recommended interval, we recommend checking out the source code of the web app or the source code of the original paper and modify it to your needs. 

In the Quantile dropdown menu, you will be presented with the option of choosing one or several quantiles. We recommend that you play with the options, but if you are looking for a comprehensive comparison of the results, choose all available quantile options.

IMPORTANT: The calculations will take some time. Increasing the size of the free parameters can significantly increase processing time. A progress bar will show, how far you are. Due to the large processing time, we recommend downloading the results, if you would like to return to it.

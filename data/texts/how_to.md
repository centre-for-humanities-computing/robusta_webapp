<br>

On Robusta you can either upload your own data for processing or use the Montecristo data used by [Herrera Malatesta and de Valeriola (2024)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0307743). <br>
The Montecristo data is meant as a test dataset, where you can play a bit around with the parameters and see potential outputs. 

**(Optional) Upload your own data** <br>
When using your own data, the dataset should be in the shapefile format and consist of a series of georeferenced points, in addition to a georeferenced polygon of the research area. <br>
Choose the tab "Upload" and upload the following shapefiles:
1. A vector file of the points to be analysed (e.g., sites, materials). You need to upload a .shp *and* a .shx file.
2. A spatial polygon file of the research area. This also needs to consist of a .shp *and* a .shx file.


**Choose number of simulations and scenarios** <br>
After specifying data (test data or own data), you need to specify the number of Monte Carlo simulations and the number of Robustness Scenarios. <br>
On the web app we recommend the following values: 
- For Monte Carlo simulations: values between 25 and 200. 
- For the Robustness Scenarios: values between 5 and 100. 

When running locally we recommend the following values: 
- For Monte Carlo simulations: values between 1000 and 5000. 
- For the Robustness Scenarios: values between 100 and 3000. 

Please note that these values are based on the standard computational power of a laptop, and increasing the values to even larger numbers can lead the app to behave in unpredicatable ways. <br>
However, users should feel free to play around with the values to see how it changes the results, as well as modifying the source code to their needs. 

**Select Quantiles** <br>
As the final step, you need to select which quantiles to run with. In the Quantile dropdown menu, you will be presented with the option of choosing one or several quantiles. 
We recommend that you play with the options, but if you are looking for a comprehensive comparison of the results, choose all available quantile options.

**IMPORTANT** <br>
Increasing the size of the free parameters can significantly increase processing time, meaning the calculations can take anything from a few minutes to several hours. 
Please be patient and do NOT refresh the page, as all data and processing will be lost. A progress bar will pop up to show how far the computations are. 
Due to the large processing time, we recommend downloading the results if you would like to return to it.

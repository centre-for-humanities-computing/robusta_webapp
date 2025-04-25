# Robusta

Robusta is an app implementing the ideas of Herrera Malatesta and de Valeriola (2024). It provides a way to quantify uncertainties in spatial analysis within the field of archaeology. 
By providing tools for data processing, visualization, and statistical analysis, Robusta enables archaeologists to make informed decisions based on spatial data while accounting for uncertainties. 
[Check out the web app for a more thorough description](https://robusta.au.dk).

---
## System Requirements and Installation

### Required Software
To use Robusta, you need the following software installed on your system:

1. **R** (Version 4.3.1 or later)
   - [Download R](https://cran.r-project.org/)
   
2. (Optional but recommended) **RStudio** 
   - [Download RStudio](https://posit.co/download/rstudio-desktop/)

### Installation Steps
Robusta is an R shiny app that opens a separate app window on your machine, which is where you can interact with the setting and the results.
To install and open the app, follow these steps:

1. Clone the repository:
In your terminal, navigate to a folder of your choice using cd and then run the following command: 
   ```bash
   git clone https://github.com/centre-for-humanities-computing/robusta_webapp.git
   ```

2. In RStudio or your preferred R environment open the Robusta Webapp project.

3. Install required dependencies 
In your R console, run the following two commands:
   ```R
   x <- c("shiny", "grid","sf", "gridExtra", "bslib", "zip", "shinyjs", "markdown",
          "spatstat", "terra", "doParallel", "foreach", "tidyverse", "scales", "ggspatial",
          "ggridges")
   
   install.packages(x)
   ```

4. Run the application by running the main script
You can either start the app by running the following command in your R console:
   ```R
   runApp("app.R")
   ```
Or you can simply pressing "Run App" in the top right corner. 

---

## How to use
A thorough guide can be seen in the app itself, but here follows a TL;DR version. 

You can either upload your own data for processing or use the Montecristi test data. 
If you change between the test data and uploading your own, we recommend refreshing the page / closing the app and re-running. 
Otherwise the app might behave in unpredictable ways when encountering an error.

**(Optional) Upload your own data** <br>
When using your own data, the dataset should be in the shapefile format and consist of a series of georeferenced points, in addition to a georeferenced polygon of the research area. <br>
Choose the tab "Upload" to upload the following shapefiles:
1. A vector file of the points to be analysed (e.g., sites, materials). You need to upload a .shp *and* a .shx file.
2. A spatial polygon file of the research area. This also needs to consist of a .shp *and* a .shx file.


**Choose number of simulations and scenarios** <br>
After specifying data, set the number of Monte Carlo simulations and the number of Robustness Scenarios. <br>
When running locally we recommend the following values: 
- For Monte Carlo simulations: values between 1000 and 5000. 
- For the Robustness Scenarios: values between 100 and 3000. 
**IMPORTANT** Dependent on the computational power of your machine, the processing can take anything from a few minutes to several hours. 

**Select Quantiles** <br>
As the final step, you need to select which quantiles to run with. In the Quantile dropdown menu, you will be presented with the option of choosing one or several quantiles. 


## Contact Information
For inquiries, support, or feedback, please reach out to:

- Eduardo Herrera Malatesta
  - Email: ehmalatesta@yahoo.com
- Sébastien de Valeriola
  - Email: sebastien.de.valeriola@ulb.be

---

## License ##
This software is [MIT licensed](./LICENSE.txt).


---

Thank you for using Robusta! Together, let's make spatial analysis in archaeology more accurate and insightful.


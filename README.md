# Robusta

Robusta is an app implementing the ideas of Herrera Malatesta and de Valeriola (2024). It provides a way to quantify uncertainties in spatial analysis within the field of archaeology. By providing tools for data processing, visualization, and statistical analysis, Robusta enables archaeologists to make informed decisions based on spatial data while accounting for uncertainties. Check out the app for a more thorough description.

---

## System Requirements and Installation

### Required Software
To use Robusta, you need the following software installed on your system:

1. **R** (Version 4.3.1 or later)
   - [Download R](https://cran.r-project.org/)
   
2. (Optional but recommended) **RStudio** 
   - [Download RStudio](https://posit.co/download/rstudio-desktop/)

### Installation Steps
1. Clone the repository:
   ```bash
   git clone https://github.com/centre-for-humanities-computing/robusta_webapp.git
   ```

2. Open the project in RStudio or your preferred R environment.

3. Install required dependencies by running the following in your R console:
   ```R
   x <- c("shiny", "grid","sf", "gridExtra", "bslib", "zip", "shinyjs", "markdown")
   
   install.packages(x)
   ```
 
4. Run the application by running the main script either via the console:
   ```R
   runApp("app.R")
   ```
   Or by simply pressing "Run App" in the top right corner. 

---

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


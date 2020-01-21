# DATR
---
Data mining with R [Test the App here!](https://kress.shinyapps.io/datr/)

---
This app will visualize Array Expression Datasets by using the [shiny framework](https://shiny.rstudio.com).

## Installation
#### There are two ways of useing the DATR app:
---
##### 1. Useing the online version
Therefore you just need to open [https://kress.shinyapps.io/datr/](https://kress.shinyapps.io/datr/) in your browser.  
The only problem here is that my computing time on the shinyapps.io server is limited to 25 hours per month, so it might happen that this page won't be reachable.

##### 2. Downloading the source code and run the app on a local device
This can be done by clicking on "Clone or download" and than on ["Download ZIP"](https://github.com/LKress/DATR/archive/master.zip).
<img src="Images/download.png" alt="drawing" width="100%"/>
After downloading the code, you should check if you have all packages installed that are needed to run the app.
The packages are:

[shiny](https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/)

[shinydashboard](https://rstudio.github.io/shinydashboard/get_started.html)

[BiocManager](https://www.bioconductor.org/install/)

[golubEsets](http://bioconductor.org/packages/release/data/experiment/html/golubEsets.html)

[RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

[DT](https://www.rdocumentation.org/packages/DT/versions/0.11)



Run this code in your R console:

```> is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])```

```> is.installed(c("shiny", "shinydashboard", "BiocManager", "golubEsets", "RColorBrewer", "DT"))```

The output should look like this:

```[1] TRUE TRUE TRUE TRUE TRUE```

If it doesn't look like this you should install the corresponding package.
I.e. if the second value isnt TRUE but FALSE you have to install shinydashboard.

If the desired output is reached change to your working directory to the directory where the downloaded source code is.

Then run this code:

```> library(shiny)```

```> runApp("DATR_App")```
 
The output should look like this:

```Listening on http://127.0.0.1:7690```

You can now open the localhost adress in your browser and you will see the app.


## The Data Tab
<img src="Images/dataTab.png" alt="drawing" width="60%"/>

On the start tab (data tab) you can find some information about the data that is visualized. The two infoboxes are showing the number of genes in the experiment and the number of patients that have been examined.
The abstract of the experiment is shown on the data tab too.

## The Diagram Tab
<img src="Images/diagramTab.png" alt="drawing" width="100%"/>

On the diagram tab you can find one slider, two selection boxes and a heatmap. 
With the slider the number of genes that should be shown in the heatmap can be regulated. In the background the app calculates the genes with the highest expression in the chosen range.

```xHighestEX = xLogarithmised[names(sort(apply(xLogarithmised,1,var), decreasing=TRUE)[1:input$numberOfGenes]),]```

This is the codeblock where the chosen number of genes with the highest expression are calculated.
The default value are 50 genes.

By choosing a distance measure in the first selection box the heatmap will be generated by useing the chosen distance measure. Availeable distance measures are: Euclidean, Maximum, Manhattan, Canberra, Binary and Minkowski.
Further information about the distance measure can be found [here](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/dist).

In the second selection box the clustering method can be chosen. This will also change the clustering method of the heatmap. Possible values are: Ward.D, Ward.D2, Single, Complete, Average (= UPGMA), Mcquitty (= WPGMA), Median (= WPGMC) and Centroid (= UPGMC).
Further information about the clustering method can be found [here](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust).



## The More Info Tab
<img src="Images/moreinfoTab.png" alt="drawing" width="100%"/>
The More Info tab shows a table of all genenames that appear in the experiment.
 
Next to the table some links for databases where these genenames can be searched in are shown. The link to this GitHub page is placed there as well.

---

A example can be testet live at [https://kress.shinyapps.io/datr/](https://kress.shinyapps.io/datr/).


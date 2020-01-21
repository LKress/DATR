# --- Load Golub Data ------

library(golubEsets)
# collect the data
data(Golub_Train)
# get the expression data
x = exprs(Golub_Train)
# get the total number of genes for the data tab
geneNr = dim(x)[1]
# sort the genenames for the table on the more info tab
genes = data.frame(sort(rownames(x)))
colnames(genes) = "Gene"
# get the patient number for the data tab
patNr = dim(x)[2]
# set all values to at least 1 to avoid NaNs
xWithoutLT1 = replace(x, x<1,1)
# logarithmize x 
xLogarithmised = log2(xWithoutLT1)

# it is possible to load other data but must be tested before use
#library(yeastCC)
#x = exprs(yeastCC)
#y = x [,5:22]
#x = replace(y,is.na(y),0)

# ------------------------

# RColorBrewer for better color of the heatmap
library("RColorBrewer")
library(shinydashboard)
# DT for inserting a table in tab: More Info
library(DT)

# user interface object
ui <- dashboardPage(
  dashboardHeader(title = "DATR"
  ),
  dashboardSidebar(
    # initializing the sidebar values 
    sidebarMenu(
      # icons can be found at https://fontawesome.com/icons?d=gallery
      menuItem("The Data", tabName = "data", icon = icon("fas fa-dna")),
      menuItem("Diagramme", tabName = "diagramme", icon = icon("fas fa-chart-pie")),
      menuItem("More Info", tabName = "moreinfo", icon = icon("fas fa-book"))
      
    
    )
  ),
  dashboardBody(
    # filling the tabs with content
    tabItems(
      tabItem(tabName = "data",
              fluidPage(
                # printing number of all genes in the experiment
                infoBox("Genes in dataset", geneNr, icon = icon("fas fa-list-ol"), fill = TRUE),
                # printing number of all patients in the experiment
                infoBox("Number of patients", patNr, icon = icon("fas fa-procedures"), fill = TRUE, color = "blue"),
                fluidRow(
                  column(12,
                    box(
                      h2("Where does the data come from?")
                    )
                  )
                ),
                fluidRow(
                  column(12,
                    box(
                      p("Golub Esets are expression data from Todd Golubs leukemia data. The dataset can be downloaded via bioconductor.", style = "font-size: 15px"),
                      a("Link to GolubEsets", href="http://bioconductor.org/packages/release/data/experiment/html/golubEsets.html"),
                      h3("Abstract of the Experiment:"),
                      # loading the abstract directly from the data so the text hasnt to be in the code and printing it in italic style
                      em(abstract(Golub_Train), style = "font-size: 16px")
                    )
                  )
                )
              )
      ),
      tabItem(tabName = "diagramme",
              sidebarLayout(
                sidebarPanel(
                  helpText("All parameters can be adjusted here. The patient IDs are on the right side of the heatmap, the gene names are at the bottom."),
                  # the user is changing the number of genes that should be shown in the heatmap
                  sliderInput(inputId = "numberOfGenes",
                              h3("Choose a number of Genes that should be shown"),
                              min = 2,
                              max = 100,
                              value = 50
                  ),
                  helpText(),
                  # the user is choosing a distance measure
                  selectInput(inputId = "distMea", 
                              h3("Distance measure"), 
                              choices = list("Euclidean" = "euclidean", "Maximum" = "maximum",
                                             "Manhattan" = "manhattan", "Canberra" = "canberra",
                                             "Binary" = "binary",  "Minkowski" = "minkowski"), selected = "euclidean"
                  ),
                  # the user is choosing the clustering method here
                  selectInput(inputId = "clustMeth", 
                              h3("Clustering method"), 
                              choices = list("Ward.D" = "ward.D", "Ward.D2" = "ward.D2",
                                             "Single" = "single", "Complete" = "complete",
                                             "Average (UPGMA)" = "average",  "Mcquitty (WPGMA)" = "mcquitty",
                                             "Median (WPGMC)" = "median",  "Centroid (UPGMC)" = "centroid"), selected = "average"
                  )
                ),
                mainPanel(
                  box(
                    h1("Heatmap Of Patients And Genes"),
                    align = "center",
                    width = 12,
                    status = "primary"
                  ),
                  box(
                    # here the heatmap will be plotted
                    plotOutput("heatmap", height = 1000),
                    width = 12
                  )
                )
              )
      ),
      tabItem(tabName = "moreinfo",
              fluidRow(
                column( 6,
                  box(width = 12,
                    h2("All Genenames to look them up"),
                    align = "center",
                    status = "primary"
                  )
                ),
                column(6,
                      box(width = 12,
                          h2("Some Links"),
                          align = "center"
                      )
                )
                
              ),
              fluidRow(
                column( width = 6,
                  # here all genes of the experiment are shown in a table for looking them up 
                  # they are sortet alphabetically
                  box( width = 12,
                    DTOutput('geneTable')
                  )
                ),
                column(6,
                       # here are some links for the user for better understanding of the data and the genes
                       box(width = 12,
                         h4("Look up the gennames at:"),
                         a("NCBI", href="https://www.ncbi.nlm.nih.gov/search/", target = "_blank"),
                         p("    "),
                         a("EBI", href="https://www.ebi.ac.uk/ebisearch/", target = "_blank"),
                         p("    "),
                         a("ENSEMBL", href="www.ensembl.org", target = "_blank"),
                         p("    "),
                         a("ENA", href="https://www.ebi.ac.uk/ena", target = "_blank"),
                         
                         h4("The code to this app:"),
                         a("GitHub", href="https://github.com/LKress/DATR", target = "_blank"),
                       ))
              )
      )
    )
  ),
  skin = c("blue")
)

# server logic unit
server <- function(input, output) {
  # rendering the heatmap plot
  output$heatmap <- renderPlot({
    # first the user chosen number of genes with the highest expression are selected
    xHighestEX = xLogarithmised[names(sort(apply(xLogarithmised,1,var), decreasing=TRUE)[1:input$numberOfGenes]),]
    # this matrix has now to be tranposed for better understanding of the heatmap
    tx = t(xHighestEX)
    # the heatmap is printet with the matrix the chosen dist and clust function and the blue color of RColorBrewer
    heatmap(tx,distfun=function(c){dist(c,method=input$distMea)}, hclustfun=function(c){hclust(c,method=input$clustMeth)}, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
  })
  # here the table on the more info tab is rendered
  output$geneTable <- renderDT(genes)
}

shinyApp(ui, server)

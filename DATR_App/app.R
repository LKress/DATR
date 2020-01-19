# --- Load Golub Data ------

library(yeastCC)
x = exprs(yeastCC)
y = x [,5:22]
x = replace(y,is.na(y),0)
#library(golubEsets)
#data(Golub_Train)
#x = exprs(Golub_Train)
geneNr = dim(x)[1]
patNr = dim(x)[2]
#x = replace(x, x<1,1)
#x = log2(x)

# ------------------------

library("RColorBrewer")
library(shinydashboard)

# user interface object
ui <- dashboardPage(
  dashboardHeader(title = "DATR"
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("DATR", tabName = "datr", icon = icon("fas fa-dna")),
      menuItem("Diagramme", tabName = "diagramme", icon = icon("fas fa-chart-pie")),
      menuItem("See also", tabName = "seealso", icon = icon("fas fa-book"))
      
    
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "diagramme",
              sidebarLayout(
                sidebarPanel(
                  helpText("All parameters can be adjusted here. The patient IDs are on the right side of the heatmap, the gene names are at the bottom."),
                  sliderInput(inputId = "numberOfGenes",
                              h3("Choose a number of Genes that should be shown"),
                              min = 2,
                              max = 100,
                              value = 50
                  ),
                  helpText(),
                  selectInput(inputId = "distMas", 
                              h3("Distance measure"), 
                              choices = list("Euclidean" = "euclidean", "Maximum" = "maximum",
                                             "Manhattan" = "manhattan", "Canberra" = "canberra",
                                             "Binary" = "binary",  "Minkowski" = "minkowski"), selected = "euclidean"
                  ),
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
                    plotOutput("heatmap", height = 1000),
                    width = 12
                  )
                )
              )
      ),
      tabItem(tabName = "datr",
              fluidPage(
                infoBox("Genes in dataset", geneNr, icon = icon("fas fa-list-ol"), fill = TRUE),
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
                      h3("Abstract of the Experiment:"),
                      em(abstract(Golub_Train), style = "font-size: 16px")
                      #em("\"Although cancer classification has improved over the past 30 years, there has been no general approach for identifying new cancer classes (class discovery) or for assigning tumors to known classes (class prediction). Here, a generic approach to cancer classification based on gene expression monitoring by DNA microarrays is described and applied to human acute leukemias as a test case. A class discovery procedure automatically discovered the distinction between acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL) without previous knowledge of these classes. An automatically derived class predictor was able to determine the class of new leukemia cases. The results demonstrate the feasibility of cancer classification based solely on gene expression monitoring and suggest a general strategy for discovering and predicting cancer classes for other types of cancer, independent of previous biological knowledge.\"", style = "font-size: 16px")
                    ),
                    box(
                      p("Experiment data"),
                      br(),
                      p("Experimenter name: Golub TR et al")    
                    )
                  )
                )
              )
      ),
      tabItem(tabName = "seealso"
      )
    )
  ),
  skin = c("blue")
)

# server logic unit
server <- function(input, output) {
  output$heatmap <- renderPlot({
    x = x[names(sort(apply(x,1,var), decreasing=TRUE)[1:input$numberOfGenes]),]
    tx = t(x)
    heatmap(tx,distfun=function(c){dist(c,method=input$distMas)}, hclustfun=function(c){hclust(c,method=input$clustMeth)}, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
    #heatmap(tx, col= colorRampPalette(brewer.pal(8, "Blues"))(25), main = "Heatmap der Patienten und Genen")
  })
}

shinyApp(ui, server)

library(shiny)
library(shinyMatrix)
library(d3r)
library(parcoords)
library(shinyBS)

# Default User Inputs ----
# Number of default traits
ntraits = 5
# Genetic variance
Gi = 0.1*diag(ntraits) + 0.9 
# Error variance
Ei = diag(ntraits) 
Ei[ntraits,ntraits] = 2
# Economic weights
#w = c(1,2) 
wi <- matrix(data = 1:ntraits, nrow = 1, ncol = ntraits)


# Define UI for application that draws in different tabs:
# 1. 3 scatterplots for up to 2 traits
# 2. 3 parallel coordinates plots for more than 2 traits
#
ui <- fluidPage(title = "Selection Indices Comparison",
                
                titlePanel("Selection Indices Comparison"),            
                
                sidebarLayout(
                  
                  sidebarPanel(

                    # Set number of traits and the size of matrices will change dynamically. 
                    sliderInput("control_num",
                                "Number of traits:",
                                min = 2, max = 10, value = ntraits),
                    # Add tooltip with instructions/info
                    bsTooltip("control_num", "Use the slider to increase or decrease the number of traits to compare",
                              "right", "hover", NULL),
                    
                    # Set population size 
                    numericInput("npop", "Set the size of the population:",
                                 min = 10, max = 1000, value = 100, step = 1), 
                    
                    # Set selection size
                    numericInput("nsel", "Set the number of selected parents:",
                                 min = 2, max = 100, value = 10, step = 1), 
                    
                    # Shiny Matrices:
                    tags$h4("Genetic variance"),
                    matrixInput(
                      "sampleG",
                      value = Gi,  
                      rows = list(
                        extend = FALSE
                      ),
                      cols = list(
                        names = FALSE
                      ),
                      class = "numeric"
                    ),
                    tags$h4("Error variance"),
                    matrixInput(
                      "sampleE",
                      value = Ei,
                      rows = list(
                        extend = FALSE
                      ),
                      cols = list(
                        names = FALSE
                      ),
                      class = "numeric"
                    ),
                    tags$h4("Economic weights"),
                    matrixInput(
                      "sampleW",
                      value = wi,
                      rows = list(
                        extend = FALSE
                      ),
                      cols = list(
                        names = FALSE
                      ),
                      class = "numeric"
                    )
                    
                  ), # endof SidebarPanel
                  # Show plots and charts 
                  mainPanel(
                    tabsetPanel(
                      tabPanel(title = "Scatterplots",
                               plotOutput("scatPlot1"),
                               plotOutput("scatPlot2"),
                               plotOutput("scatPlot3")
                      ), # endof tabPanel   
                      tabPanel(title = "Parallel Coordinates Plots",
                               parcoordsOutput("pcPlot")
                      ), # endof tabPanel
                      tabPanel(title = "Help",
                               tags$h2("Help"),
                               tags$br(),
                               tags$p("Comparing selection indices is important in breeding.
                     This application aims to help breeders explore the effect
                     of three indices in the selection of individuals."),
                               tags$p("These indices are:"),
                               tags$ol(
                                 tags$li("Base index"),
                                 tags$li("Heritability index"),
                                 tags$li("Smith-Hazel index")
                               ),
                               tags$br(),
                               tags$p("To use the app, edit the matrices on the left sidebar
                     and view how your changes affect the selection of the indices.")
                      ) # endof tabPanel             
                    )
                  ), # endof mainPanel
                ), # endof sidebarLayout
                
)


# Define server logic required to draw charts
server <- function(input, output, clientData, session) {
  
  observe({
      c_num <- input$control_num
      #updateNumericInput(session, "i_traits", value = c_num)
      updateMatrixInput(session, "sampleG", value = 0.1*diag(c_num) + 0.9)
      updateMatrixInput(session, "sampleE", value = diag(c_num))
      updateMatrixInput(session, "sampleW", value = matrix(data = 1:c_num, nrow = 1, ncol = c_num))
  })


  # Data passing function ----
  indices = function(){
    G = input$sampleG
    E = input$sampleE
    w = input$sampleW
    npop = input$npop
    nsel = input$nsel
    
    # Simulation ----
    
    # Phenotypic covariance matrix
    P = G + E
    
    # Genetic values
    g = matrix(rnorm(npop*ncol(G)),ncol=ncol(G))%*%chol(cov2cor(G))
    g = sweep(g,2,sqrt(diag(G)),"*")
    
    # Environmental (error) values
    e = matrix(rnorm(npop*ncol(E)),ncol=ncol(E))%*%chol(cov2cor(E))
    e = sweep(e,2,sqrt(diag(E)),"*")
    
    # Phenotype
    p = g+e
    
    # Select using base index
    take_base = order(p%*%w[1,], decreasing=TRUE)[1:nsel]
    p_base = p[take_base,]
    g_base = g[take_base,]
    
    # Select using h2 index
    h = w[1,] * diag(G) / diag(P)
    take_h2 = order(p%*%h, decreasing=TRUE)[1:nsel]
    p_h2 = p[take_h2,]
    g_h2 = g[take_h2,]
    
    # Select using Smith-Hazel index
    b = solve(P) %*% G %*% w[1,]
    take_smith = order(p%*%b, decreasing=TRUE)[1:nsel]
    p_smith = p[take_smith,]
    g_smith = g[take_smith,]
    
    # Make plots ----
    # Only works for 2 traits
    
    # Mean g
    g_mean = colMeans(g)
    g_base_mean = colMeans(g_base)
    g_h2_mean = colMeans(g_h2)
    g_smith_mean = colMeans(g_smith)
    
    # Index means
    i_init = mean(g%*%w[1,])
    i_base = mean(g_base%*%w[1,]) - i_init
    i_h2 = mean(g_h2%*%w[1,]) - i_init
    i_smith = mean(g_smith%*%w[1,]) - i_init
    
    # Extreme plot point
    lim = max(abs(p))
    
    return(list
           (
             "p"=p,
             "take_base" = take_base,
             "take_h2" = take_h2,
             "take_smith" = take_smith,
             "p_base"=p_base,
             "p_h2"=p_h2,
             "p_smith"=p_smith,
             "i_base"=i_base,
             "i_h2"=i_h2,
             "i_smith"=i_smith,
             "lim"=lim,
             "g_mean"=g_mean,
             "g_base_mean"=g_base_mean,
             "g_h2_mean"=g_h2_mean,
             "g_smith_mean"=g_smith_mean,
             "npop"=npop,
             "nsel"=nsel
             )
           )
  } 

  
  # First Tab 
  output$scatPlot1 <- renderPlot({

    data <- indices()
      p = data$p
      take_base = data$take_base
      take_h2 = data$take_h2
      take_smith = data$take_smith
      p_base = data$p_base
      p_h2 = data$p_h2
      p_smith = data$p_smith
      i_base = data$i_base
      i_h2 = data$i_h2
      i_smith = data$i_smith
      lim = data$lim
      g_mean = data$g_mean
      g_base_mean = data$g_base_mean
      g_h2_mean=data$g_h2_mean
      g_smith_mean=data$g_smith_mean 

    
    op = par(mfrow=c(1,3))
    # Base index
    plot(p,
         xlab="Trait 1",
         ylab="Trait 2",
         main=paste0("Base Index (gain=",
                     round(i_base, 2),")"),
         pch=20,
         col="lightblue",
         xlim=c(-lim,lim),
         ylim=c(-lim,lim))
    points(p_base,
           pch=20,col="pink")
    points(g_mean[1],g_mean[2],
           pch=15,cex=1.5)
    points(g_base_mean[1],
           g_base_mean[2],
           pch=17,cex=2,col="red")
    
    
    
    # Heritability index
    plot(p,
         xlab="Trait 1",
         ylab="Trait 2",
         main=paste0("Heritability Index (gain=",
                     round(i_h2, 2),")"),
         pch=20,
         col="lightblue",
         xlim=c(-lim,lim),
         ylim=c(-lim,lim))
    points(p_h2,
           pch=20,col="pink")
    points(g_mean[1],g_mean[2],
           pch=15,cex=1.5)
    points(g_h2_mean[1],
           g_h2_mean[2],
           pch=17,cex=2,col="red")
    
    # Smith-Hazel index
    plot(p,
         xlab="Trait 1",
         ylab="Trait 2",
         main=paste0("Smith-Hazel Index (gain=",
                     round(i_smith, 2),")"),
         pch=20,
         col="lightblue",
         xlim=c(-lim,lim),
         ylim=c(-lim,lim))
    points(p_smith,
           pch=20,col="pink")
    points(g_mean[1],g_mean[2],
           pch=15,cex=1.5)
    points(g_smith_mean[1],
           g_smith_mean[2],
           pch=17,cex=2,col="red")
    
    par(op)
    
    
    
  })    
  
  # Second Tab    
  output$pcPlot <- renderParcoords({
    
    data <- indices()
      p = data$p
      take_base = data$take_base
      take_h2 = data$take_h2
      take_smith = data$take_smith
      p_base = data$p_base
      p_h2 = data$p_h2
      p_smith = data$p_smith
      i_base = data$i_base
      i_h2 = data$i_h2
      i_smith = data$i_smith
      lim = data$lim
      g_mean = data$g_mean
      g_base_mean = data$g_base_mean
      g_h2_mean=data$g_h2_mean
      g_smith_mean=data$g_smith_mean 
      npop=data$npop
      nsel=data$nsel
   
    #TV: data manipulation to get not selected p
    take_all<-c(take_base,take_h2,take_smith)
    take_all<-unique(take_all)
    hund<-c(1:npop)
    hund<-hund[ ! hund %in% take_all ]
    pns=p[hund,]
    # TV TODO : add -+lim entries to normalise axes
    maxl <- g_mean
    for (i in 1:length(maxl)) { maxl[i] <- lim}
    minl <- g_mean
    for (i in 1:length(minl)) { minl[i] <- -lim}
    mm<-rbind(minl,maxl)
    rownames(mm)<-c(npop+1,npop+2)
    pns<-rbind(pns,mm)
    #TV: Transform matrix to data frame
    pns <- as.data.frame(pns) 
    #TV: Add column $selection "not_selected"
    pns$selection <- "not selected"
    #TV: Do the same for selected
    g_mean <- as.data.frame(t(g_mean))
    g_mean$selection <- "Mean"
    p_base <- as.data.frame(p_base)
    p_base$selection <- "Base"
    g_base_mean <- as.data.frame(t(g_base_mean))
    g_base_mean$selection <- "Base Mean"
    p_h2 <- as.data.frame(p_h2)
    p_h2$selection <- "H2"
    g_h2_mean <- as.data.frame(t(g_h2_mean))
    g_h2_mean$selection <- "H2 Mean"
    p_smith <- as.data.frame(p_smith)
    p_smith$selection <- "Smith"
    g_smith_mean <- as.data.frame(t(g_smith_mean))
    g_smith_mean$selection <- "Smith Mean"
    
    pd<-rbind(pns,g_mean,p_base,g_base_mean,p_h2,g_h2_mean,p_smith,g_smith_mean)
    row.names(pd) <- NULL

    #p<-rbind(pns,p_base,p_h2,p_smith)
    #op = par(mfrow=c(1,3))
    


    parcoords(
      pd,
      brushMode = '1D-axes-multi', # "1D-axes", "1D-axes-multi", or "2D-strums" 
      reorderable = TRUE,
      color = list(
         # discrete or categorical column
         colorScale = "scaleOrdinal",
         colorBy = "selection",
         colorScheme = "schemePaired" # "schemeDark2" # "schemeCategory10"
       ),
      autoresize = TRUE,
      withD3 = TRUE,
      alphaOnBrushed = 0, # 0.1,
      height = 400
    )
    
   
    
  })   
  

}

# Run the application 
shinyApp(ui = ui, server = server)

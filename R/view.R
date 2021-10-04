#' Visulization of rMAUPS results
#'
#' @docType methods
#' @name view
#' @rdname view
#'
#' @param outdir The path to rMAUPS results.
#'
#' @author Wubing Zhang
#'
#' @return Open the shiny app.
#' @import shiny
#' @export
#'
view <- function(outdir = "./"){
  requireNamespace("shiny")
  ui <- fluidPage(
    headerPanel('Visualization of rMAUPS results'),
    textInput(inputId = "outdir", label = "Enter the output directory",
              value = outdir, width = 10000),
    actionButton(inputId = "submit", label = "submit"),
    h1(),
    tabsetPanel(
      tabPanel("Quality control", h4("Three type of quality controls:
               \n\t1. The distribution of protein abundance for each sample.
               \n\t2. Principal component analysis \n\t3. Pairwise correlation of proteins within the same complexes
               \n\t4. The distribution of dropouts (when there are NAs in the data)."),
               uiOutput("qc")),
      tabPanel("Imputation - QC", h4("After imputation using knn, the quality control analyses are done again."),
               "*hint: No output when there is no missing value in the data.",
               uiOutput("impute")),
      tabPanel("Comparison", h4("Differential expression analysis is performed for each comparison at three different levels,
                including protein level, pathway level, and protein complex level."), uiOutput("compare")),
      tabPanel("Integration", h4("The same comparison in different datasets are combined and visualized in this section."),
               uiOutput("integrate"))
    )
  )

  server <- function(input, output){
    pnglist <- eventReactive(input$submit, {
      readRDS(paste0(input$outdir,"/", "summary_list.rds")) })

    output$qc <- renderUI({
      plot_output_list <- lapply(1:length(pnglist()$qc), function(i) {
        plotname <- paste("qc", i, sep="")
        plotOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, c(plot_output_list))
    })
    output$impute <- renderUI({
      plot_output_list <- lapply(1:length(pnglist()$impute), function(i) {
        plotname <- paste("impute", i, sep="")
        plotOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, c(plot_output_list))
    })
    output$compare <- renderUI({
      plot_output_list <- lapply(1:length(pnglist()$compare), function(i) {
        plotname <- paste("compare", i, sep="")
        plotOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, c(plot_output_list))
    })
    output$integrate <- renderUI({
      plot_output_list <- lapply(1:length(pnglist()$integrate), function(i) {
        plotname <- paste("integrate", i, sep="")
        plotOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, c(plot_output_list))
    })

    #### renderplot all the figures
    for (i in 1:20) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i ## Very necessary
        output[[paste("qc", my_i, sep="")]] <- renderPlot({ pnglist()$qc[[my_i]]})
        output[[paste("impute", my_i, sep="")]] <- renderPlot({ pnglist()$impute[[my_i]]})
        output[[paste("compare", my_i, sep="")]] <- renderPlot({ pnglist()$compare[[my_i]]})
        output[[paste("integrate", my_i, sep="")]] <- renderPlot({ pnglist()$integrate[[my_i]]})
      })
    }
  }
  shinyApp(ui = ui, server = server)
}


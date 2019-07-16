library(shiny)
library(shinyMatrix)

ui <- fluidPage(
  pageWithSidebar(
    headerPanel('spBayes Conjugate Linear Model'),
    sidebarPanel(
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File. Make sure the CSV file selected has the response variable outcomes listed in the first column
              followed by columns corresponding to each of the explanatory variable values.",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),

      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),

      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),


      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),

      tags$hr(),

      selectInput(inputId = "select", "Select variables to regress on", choices = " ", multiple = TRUE),
      numericInput(inputId = "nsamp", 'Number of Samples', value = 10),
      textInput(inputId = "mu_beta", "beta.prior.mean: Enter a vector (comma delimited)", value = c()),
      numericInput(inputId = "prior_shape", "Prior Shape", value = -3.5),
      numericInput(inputId = "prior_rate", "Prior Rate", value = 0),

      tags$hr(),


      actionButton("go", "Submit")
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Sigma Squared Samples", plotOutput('hist1')),
        #tabPanel("Beta Samples", uiOutput("plots")),
        tabPanel("Tables", textOutput('txt'), dataTableOutput('table1'), dataTableOutput('table2'))
      )
    )
  )
)

server <- function(session, input, output){


  data <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = input$header, sep = input$sep)
  })

  filtereddata <- eventReactive({
    input$update
    data()},  {
    req(data())
    if(is.null(input$select) || input$select == "")
      data() else
        {data[, colnames(data()) %in% input$select] }

  })


  observeEvent(data(), {
    updateSelectInput(session, "select", choices=colnames(data()[,2:ncol(data())]))
  })


  observeEvent(input$select, {
  req(input$select)
    vec = rep(0, times = length(input$select)+1)
    updateTextInput(session, "mu_beta", value = vec)
  })


  Y <- reactive({names(data())[1]})

  # beta.prior.mean <- reactive({rep(0, times = (length(input$select) + 1))}) # mu_beta
  beta.prior.precision <- reactive({matrix(0, nrow=(length(input$select) + 1),
                                           ncol= (length(input$select) + 1))}) # V_beta


  m.1 <- eventReactive(input$go, {
    req(input$mu_beta)
    req(input$select)
    spBayes::bayesLMConjugate(formula = as.formula(paste(Y(), "~", paste(unlist(strsplit(input$select, " ")), collapse = "+"))),
                              data = data(),
                              n.samples = input$nsamp, beta.prior.mean = as.numeric(unlist(strsplit(input$mu_beta, ","))),
                              beta.prior.precision = beta.prior.precision(),
                              prior.shape = input$prior_shape, prior.rate = input$prior_rate)
  })


  output$txt <- renderText("The first table displays the empirical means, the standard deviations,
                            and the standard errors for each variable.
                           The second table displays the quantiles for each variable.")
  output$table1 <- renderDataTable(round(summary(m.1()$p.beta.tauSq.samples)[[1]], 3),
                                   options = list(include.rownames = TRUE, initComplete = JS(
                                     "function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                     "}")))
  output$table2 <- renderDataTable(round(summary(m.1()$p.beta.tauSq.samples)[[2]], 3),
                                   options = list(include.rownames = TRUE, initComplete = JS(
                                     "function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                     "}")))

  output$hist1 <- renderPlot({
    sigma_sq <- m.1()$p.beta.tauSq.samples[,ncol(m.1()$p.beta.tauSq.samples)]
    hist(sigma_sq, main = "Sigma Squared Samples")}
  )

  # output$hist2 <- renderPlot({
  #   # beta_vals <- m.1()$p.beta.tauSq.samples[,2]
  #   # hist(beta_vals, main = "Beta Samples")}
  #   for(i in 1:(ncol(m.1()$p.beta.tauSq.samples)-1)){
  #     beta_vals <- m.1()$p.beta.tauSq.samples[,i]
  #     hist(beta_vals, main = paste0("Beta", i, " Samples"))
  #   }
  # })

  # output$plots <- renderUI({
  #   plot_output_list <- lapply(1:(ncol(m.1()$p.beta.tauSq.samples)-1), function(i) {
  #     plotname <- paste0("plot", i)
  #     htmlOutput(plotname)
  #   })
  #
  #   # Convert the list to a tagList - this is necessary for the list of items
  #   # to display properly.
  #   tagList(plot_output_list)
  # })
  #
  # eventReactive(data(), {
  #   req(input$file1)
  #   for (i in 1:(ncol(m.1()$p.beta.tauSq.samples)-1)) {
  #   local({
  #     my_i <- i
  #     plotname <- paste("Beta", my_i, " Samples")
  #
  #     output[[plotname]] <- renderPlot({
  #       beta_vals <- m.1()$p.beta.tauSq.samples[,my_i]
  #       hist(beta_vals(), main = plotname)
  #     })
  #   })
  # }})



}

shinyApp(ui, server)

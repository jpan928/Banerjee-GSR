library(shiny)
library(spBayes)
library(DT)
data("FORMGMT.dat")

ui <- fluidPage(
  pageWithSidebar(
    headerPanel('spBayes Geostat Exact using Covariates'),
    sidebarPanel(

      # Input: Select a file ----
      fileInput("file1", "Choose CSV File. Make sure the CSV file selected
                has the response variable outcomes listed in the first column
                followed by columns corresponding to each of the explanatory
                variable values.",
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

      numericInput(inputId = "nsamp", 'Number of Samples', value = 10),
      # checkboxGroupInput(inputId = 'xvar', 'Explanatory Variables',
      #                    choices = names(FORMGMT.dat)[2:7], selected = names(FORMGMT.dat)[2:5]),
      # checkboxGroupInput(inputId = 'xvar', 'Explanatory Variables',
      #                    choices = "", selected = ""),
      selectInput(inputId = "xvar", "Explanatory Variables", choices = " ", multiple = TRUE),
      textInput(inputId = "mu_beta", "beta.prior.mean: Enter a vector (comma delimited)", value = c()),
      # selectInput(inputId = 'yvar', 'Response Variable', names(FORMGMT.dat)[1]),
      numericInput(inputId = "phi", "phi", value = 0.0012),
      numericInput(inputId = "prior_shape", "Prior Shape", value = -3.5),
      numericInput(inputId = "prior_rate", "Prior Rate", value = 0),
      numericInput(inputId = "alpha", "alpha", value = 1/1.5), # what is this used for?
      numericInput(inputId = "prior_shape", "Sigma Squared Prior Shape", value = 2.0),
      numericInput(inputId = "prior_rate", "Sigma Squared Prior Rate", value = 10.0),

      tags$hr(),

      actionButton("go", "Submit")
    ),
    mainPanel(
      tabsetPanel( type = "tabs",
                   tabPanel("Plots", plotOutput('hist1'),
                            plotOutput('hist2')),
                   tabPanel("Table", textOutput('txt'), dataTableOutput('table')))
    )
  )
)



server <- function(input, output, session) {

  data <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = input$header, sep = input$sep)
  })

  filtereddata <- eventReactive({
    input$update
    data()},  {
      req(data())
      if(is.null(input$xvar) || input$xvar == "")
        data() else
        {data[, colnames(data()) %in% input$xvar] }

    })

  observeEvent(data(), {
    updateSelectInput(session, "xvar", choices=colnames(data()[,grep("X.", colnames(data()))]))
  })


  Y <- eventReactive(data(),
                   {names(data()[grep("Y", colnames(data()))])})

  p <- eventReactive(input$xvar,
                    {req(input$file1)
                      length(input$xvar) + 1})

  observeEvent(p(), {
    vec = as.matrix(rep(0, times = p()))
    updateTextInput(session, "mu_beta", value = vec)
  })

  coords <- eventReactive(data(),
                          {(cbind(data()[,grep("Long", colnames(data()))],
                                  data()[,grep("Lat", colnames(data()))]))*(pi/180)*6378})

  beta.prior.precision <- reactive({matrix(0, nrow = p(), ncol = p())}) # V_beta

  m.1 <- eventReactive(input$go, {
    req(input$xvar)
    req(input$file1)
    req(input$mu_beta)
    spBayes::bayesGeostatExact(as.formula(paste(Y(), "~", paste(unlist(strsplit(input$xvar, " ")), collapse = "+"))),
                               data = data(),
                               n.samples=input$nsamp,
                               beta.prior.mean=as.numeric(unlist(strsplit(input$mu_beta, ","))),
                               beta.prior.precision=beta.prior.precision(),
                               coords=coords(), phi=input$phi, alpha=input$alpha,
                               sigma.sq.prior.shape=input$prior_shape,
                               sigma.sq.prior.rate=input$prior_rate)
  })


  output$txt <- renderText({"The following table displays the posterior samples for the defined
    parameters."})
   output$table <- renderDataTable({as.data.frame(round(m.1()$p.samples, 3))},
                                   options = list(include.rownames = TRUE,
                                                  initComplete = JS("function(settings, json) {",
                                                                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                                                    "}")))


  output$hist1 <- renderPlot({
    sigma_sq <- m.1()$p.samples[,2]
    hist(sigma_sq, main = "Sigma Squared Samples")}
  )

  output$hist2 <- renderPlot({
    tau_sq <- m.1()$p.samples[,3]
    hist(tau_sq, main = "Tau Squared Samples")}
  )

}

shinyApp(ui = ui, server = server)


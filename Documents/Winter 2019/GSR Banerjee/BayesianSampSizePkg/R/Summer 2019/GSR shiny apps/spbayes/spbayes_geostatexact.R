library(shiny)
library(spBayes)
library(MBA)
library(DT)

# Use this dataset as an example. Save it into your computer and then upload it as a CSV in this app
data("FBC07.dat")

ui <- fluidPage(
  pageWithSidebar(
    headerPanel('spBayes Geostat Exact (modeling with Intercept only)'),

    sidebarPanel(
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File. Make sure the CSV file selected follows the following format:

              1) Columns corresponding to the
              response variables are labeled as 'Y.number'

              2) Columns corresponding to the x and y coordinates are
              labeled as coord.X and coord.Y respectively.",
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

      sliderInput(inputId = "range", 'Range of Interest', value = c(1,10), min = 1, max = 100),
      selectInput(inputId = "select", "Select response variable of interest", choices = " ", multiple = FALSE),
      numericInput(inputId = "nsamp", 'Number of Samples', value = 10),
      textInput(inputId = "mu_beta", "beta.prior.mean: Enter a vector (comma delimited)", value = c()),
      numericInput(inputId = "phi", "phi", value = 0.15),
      numericInput(inputId = "nu", "nu", value = 0.5), # what is this used for?
      numericInput(inputId = "alpha", "alpha", value = 1), # what is this used for?
      numericInput(inputId = "p", "p", value = 1),
      numericInput(inputId = "prior_shape", "Sigma Squared Prior Shape", value = 2.0),
      numericInput(inputId = "prior_rate", "Sigma Squared Prior Rate", value = 5.0),

      tags$hr(),


      actionButton("go", "Submit")
    ),
    mainPanel(

      tabsetPanel( type = "tabs",
                   tabPanel("Plots", plotOutput('hist1'),
                            plotOutput('hist2')),
                   tabPanel("Tables", textOutput('txt'), dataTableOutput('table1'), dataTableOutput('table2')),
                   tabPanel("Contour Plots", plotOutput('img')))
    )
)
)


server <- function(input, output, session) {

  d <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = input$header, sep = input$sep)
  })

  data <- eventReactive(d(), {
    na.omit(d())
  })

  observeEvent(data(), {
    updateSelectInput(session, "select", choices=colnames(data()[,grep("Y.", colnames(data()))]))
  })

  observeEvent(data(),
    {updateSliderInput(session, "range", min = 1, max = nrow(data()))})

  Y <- eventReactive(input$range,
    {req(input$file1)
      req(input$select)
      data()[input$range[1]:input$range[2], paste0("", input$select, "")]})


  coords <- eventReactive(input$range,
      {req(input$file1)
      as.matrix(data()[input$range[1]:input$range[2], c("coord.X", "coord.Y")])})


  observeEvent(input$p, {
    req(input$p)
    vec = as.matrix(rep(0, times = input$p))
    updateTextInput(session, "mu_beta", value = vec)
  })


  beta.prior.precision <- reactive({matrix(0, nrow=input$p, ncol= input$p)}) # V_beta


  m.1 <- eventReactive(input$go, {
    req(input$mu_beta)
    req(input$range)
    spBayes::bayesGeostatExact(Y() ~ 1, n.samples=input$nsamp,
                               beta.prior.mean=as.numeric(unlist(strsplit(input$mu_beta, ","))),
                               beta.prior.precision=beta.prior.precision(),
                               coords=coords(), phi=input$phi, alpha=input$alpha,
                               sigma.sq.prior.shape=input$prior_shape,
                               sigma.sq.prior.rate=input$prior_rate)
  })


  obs.surf <- eventReactive(input$go, {
    par(mfrow=c(1,2))
    mba.surf(cbind(coords(), Y()), no.X=100, no.Y=100, extend=T)$xyz.est})

  w.hat <- eventReactive(input$go, {rowMeans(m.1()$sp.effects)})
  w.surf <-
    eventReactive(input$go, {mba.surf(cbind(coords(), w.hat()), no.X=100, no.Y=100, extend=T)$xyz.est})


  output$img <- renderPlot({image(obs.surf(), xaxs = "r", yaxs = "r", main="Observed response")
    points(coords())
    contour(obs.surf(), add=T)

    image(w.surf(), xaxs = "r", yaxs = "r", main="Estimated random effects")
    points(coords())
    contour(w.surf(), add=T)}
    )

  output$txt <- renderText({"The first table displays the empirical means, standard
  deviations and standard errors of the means for each variable. The second
  table displays the quantiles of the variables."})
  output$table1 <- renderDataTable({round(summary(m.1()$p.samples)[[1]], 3)},
                                   options = list(include.rownames = TRUE, initComplete = JS(
                                     "function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                     "}")))

  output$table2 <- renderDataTable({round(summary(m.1()$p.samples)[[2]], 3)},
                                  options = list(include.rownames = TRUE, initComplete = JS(
                                    "function(settings, json) {",
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


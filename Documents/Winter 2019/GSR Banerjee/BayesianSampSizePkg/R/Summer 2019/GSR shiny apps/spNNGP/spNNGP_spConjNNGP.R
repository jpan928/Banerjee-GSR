# The function spConjNNGP fits Gaussian univariate Bayesian conjugate spatial regression models using Nearest Neighbor Gaussian Processes (NNGP).
library(shiny)
library(spNNGP)


##Fit a Conjugate NNGP model and predict for the holdout
sigma.sq.IG <- c(2, sigma.sq)
g <- 10
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))
colnames(theta.alpha) <- c("phi", "alpha")


ui <- fluidPage(
  pageWithSidebar(
    headerPanel('spNNGP::spConjNNGP'),
    sidebarPanel(
      selectInput(inputId = "select", label = "Choose whether you wish to generate data or upload your own dataset",
                  choices = c("Generate data", "Upload dataset"), selected = NULL),

      tags$hr(),

      conditionalPanel(

        condition = "input.select == 'Generate data'",

        numericInput(inputId = "n", label = "Number of data observations to generate", value = 100),
        numericInput(inputId = "sig.sq", label = "Sigma Squared", value = 5),
        numericInput(inputId = "tau.sq", label = "Tau Squared", value = 1),
        numericInput(inputId = "phi", label = "Phi", value = 3/0.5),
        selectInput(inputId = 'score.rule', 'Score Rule', choices = c('crps', 'rmspe'), multiple = FALSE, selected = 'crps'),
        textInput(inputId = "B", "B matrix: Enter a vector (comma delimited)", value = "1, 5"),
        numericInput(inputId = "n.neighbors", 'Number of Neighbors', value = 10),
        numericInput(inputId = "k.fold", "Number of k-folds for cross-validation", value = 5),
        numericInput(inputId = 'n.omp.threads', 'Number of threads to use for SMP parallel processing', value = 1),
        selectInput(inputId = 'cov.model', "Covariance Model", choices = c('exponential', 'spherical', 'gaussian'), selected = 'exponential',
                    multiple = FALSE),
        numericInput(inputId = 'n.samples', 'Number of Samples', value = 500),
        numericInput(inputId = 'n.report', 'Interval to report sampling progress', value = 500),
        actionButton("go", "Submit")
      ),

      conditionalPanel(
        condition = "input.select == 'Upload dataset'",
        # Input: Select a file ----
        fileInput("file1", "Choose CSV File. Be sure there are columns corresponding to the data observations (y), the
                covariate (x), and the coordinates (with columns labeled as coords1 and coords2).",
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

        numericInput(inputId = "n.neighbors2", 'Number of Neighbors', value = 10),
        numericInput(inputId = "phi2", label = "Phi", value = 3/0.5),
        numericInput(inputId = "sig.sq2", label = "Sigma Squared", value = 5),
        numericInput(inputId = "tau.sq2", label = "Tau Squared", value = 1),
        selectInput(inputId = 'score.rule2', 'Score Rule', choices = c('crps', 'rmspe'), multiple = FALSE, selected = 'crps'),
        numericInput(inputId = "k.fold2", "Number of k-folds for cross-validation", value = 5),
        numericInput(inputId = 'n.omp.threads2', 'Number of threads to use for SMP parallel processing', value = 1),
        selectInput(inputId = 'cov.model2', "Covariance Model", choices = c('exponential', 'spherical', 'gaussian'), selected = 'exponential',
                    multiple = FALSE),
        numericInput(inputId = 'n.samples2', 'Number of Samples', value = 500),
        numericInput(inputId = 'n.report2', 'Interval to report sampling progress', value = 500),
        actionButton("go2", "Submit")
      )
    ),

    mainPanel(

     tabsetPanel(type = "tabs",
                tabPanel("Regression Coefficient Estimates", textOutput('txt'), dataTableOutput('table1')),
                tabPanel("K-Fold Cross Validation", textOutput('txt2'), dataTableOutput('table2')))
  )
)
)





server <- function(input, output, session) {


  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }


  data <- reactive({
    req(input$file1)
    req(input$select)
    if(input$select == "Upload dataset"){
      read.csv(input$file1$datapath, header = input$header, sep = input$sep)
    }
  })


  ##################
  # Upload Dataset #
  ##################
  y.up.data <- eventReactive(data(), {
    req(input$file1)
    req(input$select)
    if(input$select == "Upload dataset"){
      data()[,which(colnames(data()) == "y")]
    }
  })


  coords.up.data <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      cbind(data()[,which(colnames(data()) == "coords1")],
            data()[,which(colnames(data()) == "coords2")])
    }
  })

  x.up.data <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      cbind(1, data()[,which(colnames(data()) == "x")])
    }
  })

  D.upload <- eventReactive(data(),{
    req(input$select)
    if(input$select == "Upload dataset"){
      as.matrix(dist(coords.up.data()))
    }
  })

  R.upload <- eventReactive(data(), {
    exp(-input$phi2 * D.upload())
  })

  w.upload <- eventReactive(data(),{
    req(input$select)
    if(input$select == "Upload dataset"){
      rmvn(1, rep(0, nrow(data())), input$sig.sq2 * R.upload())
    }
  })

  ho.upload <- eventReactive(data(), {
    req(input$file1)
    req(input$select)
    #req(input$n.upload)
    if(input$select == "Upload dataset"){
      sample(1:nrow(data()), nrow(data())/2)
    }
  })

  y.ho.upload <-  eventReactive(data(), {
    req(input$file1)
    req(input$select)
    if(input$select == "Upload dataset"){
      y.up.data()[ho.upload()]
    }
  })

  x.ho.upload <-  eventReactive(data(), {
    req(input$file1)
    req(input$select)
    # req(input$n.upload)
    if(input$select == "Upload dataset"){
      x.up.data()[ho.upload(),,drop=FALSE]
    }
  })

  w.ho.upload <-  eventReactive(data(), {
    w.upload()[ho.upload()]
  })

  coords.ho.upload <-  eventReactive(data(), {
    req(input$file1)
    req(input$select)
    if(input$select == "Upload dataset"){
      coords.up.data()[ho.upload(),]
    }
  })





  ####################
  # Generate Dataset #
  ####################
  coords.gen.data <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      cbind(runif(input$n,0,1), runif(input$n,0,1))
    }
  })

  x.gen.data <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      cbind(1, rnorm(input$n))
    }
  })

  D <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      as.matrix(dist(coords.gen.data()))
    }
  })

  R <- reactive({
    exp(-input$phi * D())
  })

  w <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      rmvn(1, rep(0, input$n), input$sig.sq * R())
    }
  })

  y.gen.data <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      rnorm(input$n, x.gen.data()%*%as.matrix(as.numeric(unlist(strsplit(input$B, ",")))) + w(), sqrt(input$tau.sq))
    }
  })

  ho <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      sample(1:input$n, 50)
    }
  })



  y.ho <-  reactive({
    req(input$select)
    if(input$select == "Generate data"){
      y.gen.data()[ho()]
    }
  })

  x.ho <-  reactive({
    req(input$select)
    if(input$select == "Generate data"){
      x.gen.data()[ho(),,drop=FALSE]
    }
  })

  w.ho <-  reactive({
    w()[ho()]
  })

  coords.ho <-  reactive({
    req(input$select)
    if(input$select == "Generate data"){
      coords.gen.data()[ho(),]
    }
  })


  y.new.1 <-  eventReactive(input$select, {
    if(input$select == "Generate data"){
      y.gen.data()[-ho()]
    }
  })
  y.new.2 <-  eventReactive(data(), {
    req(input$file1)
    if(input$select == "Upload dataset"){
      y.up.data()[-ho.upload()]
    }
  })


  x.new.1 <-  eventReactive(input$select, {
    if(input$select == "Generate data"){
      x.gen.data()[-ho(),,drop = FALSE]
    }
  })
  x.new.2 <-  eventReactive(data(), {
    if(input$select == "Upload dataset"){
      req(input$file1)
      x.up.data()[-ho.upload(),,drop = FALSE]
    }
  })


  w.new.1 <- eventReactive(input$select, {
    if(input$select == "Generate data"){
      w()[-ho(),,drop = FALSE]
    }
  })
  w.new.2 <- eventReactive(data(), {
    if(input$select == "Upload dataset"){
      req(input$file1)
      w.upload()[-ho.upload(),,drop = FALSE]
    }
  })


  coords.new.1 <- eventReactive(input$select, {
    if(input$select == "Generate data"){
      coords.gen.data()[-ho(),]
    }
  })
  coords.new.2 <-  eventReactive(data(), {
    if(input$select == "Upload dataset"){
      req(input$file1)
      coords.up.data()[-ho.upload(),]
    }
  })



  m.c.1 <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spNNGP::spConjNNGP(y.new.1()~x.new.1()-1, coords=coords.new.1(), n.neighbors = input$n.neighbors,
                         X.0 = x.ho(), coords.0 = coords.ho(),
                         k.fold = input$k.fold, score.rule = input$score.rule,
                         n.omp.threads = input$n.omp.threads,
                         theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = input$cov.model)
    }
  })
  m.c.2 <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spNNGP::spConjNNGP(y.new.2()~x.new.2()-1, coords=coords.new.2(), n.neighbors = input$n.neighbors2,
                         X.0 = x.ho.upload(), coords.0 = coords.ho.upload(),
                         k.fold = input$k.fold2, score.rule = input$score.rule2,
                         n.omp.threads = input$n.omp.threads2,
                         theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = input$cov.model2)
    }
  })


  output$txt <- renderText({"The table below shows the regression coefficient estimates corresponding to the returned theta.alpha."})

  output$table1 <- renderDataTable({
    if(input$select == "Generate data"){
      round(m.c.1()$beta.hat,3)
    }else if(input$select == "Upload dataset"){
      round(m.c.2()$beta.hat,3)
    }
  }, options = list(columns = list(
    list(title = "x1"), list(title = "x2")),
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
      "}"))
  )

  output$txt2 <- renderText({"The table below shows results from the k-fold cross validation."})

  output$table2 <- renderDataTable({
    if(input$select == "Generate data"){
      round(m.c.1()$k.fold.scores,3)
    }else if(input$select == "Upload dataset"){
      round(m.c.2()$k.fold.scores,3)
    }
  }, options = list(
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
      "}")))

}

shinyApp(ui = ui, server = server)


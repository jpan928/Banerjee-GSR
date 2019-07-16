##Fit a Response and Sequential NNGP model
library(shiny)
library(spNNGP)

# Future edits could make this more customizable
phi = 3 / 0.5
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))


ui <- fluidPage(
  pageWithSidebar(
  headerPanel('spNNGP::spPredict'),
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
    plotOutput('plot1'),
    plotOutput('plot2')
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

  # starting <- eventReactive(data(), {
  #   req(input$file1)
  #   req(input$select)
  #   if(input$select == "Upload dataset"){
  #     list("phi"=input$phi2, "sigma.sq"=input$sig.sq2, "tau.sq"=input$tau.sq2)
  #   }
  # })
  #
  # tuning <- eventReactive(data(), {
  #   req(input$file1)
  #   req(input$select)
  #   if(input$select == "Upload dataset"){
  #     list("phi"=input$phi2, "sigma.sq"=input$sig.sq2, "tau.sq"=input$tau.sq2)
  #   }
  # })


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


  m.s.1 <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spNNGP(y.new.1()~x.new.1()-1, coords=coords.new.1(), starting=starting, method='sequential', n.neighbors=input$n.neighbors,
             tuning=tuning, priors=priors, cov.model=input$cov.model,
             n.samples=input$n.samples, n.omp.threads=input$n.omp.threads, n.report=input$n.report)
    }
  })

  m.s.2 <- eventReactive(input$go2, {
    req(data())
    #req(input$select)
    if(input$select == "Upload dataset"){
      spNNGP(y.new.2()~x.new.2()-1, coords=coords.new.2(), starting=starting, method='sequential', n.neighbors=input$n.neighbors2,
             tuning=tuning, priors=priors, cov.model=input$cov.model2,
             n.samples=input$n.samples2, n.omp.threads=input$n.omp.threads2, n.report=input$n.report2)
    }
  })

  p.s.1 <- eventReactive(input$go, {
    if(input$select == "Generate data"){
      spNNGP::spPredict(m.s.1(), X.0 = x.ho(), coords.0 = coords.ho(), n.omp.threads=input$n.omp.threads)
    }
  })


  p.s.2 <- eventReactive(input$go2, {
    req(data())
    #req(input$select)
    if(input$select == "Upload dataset"){
      spNNGP::spPredict(m.s.2(), X.0 = x.ho.upload(), coords.0 = coords.ho.upload(), n.omp.threads=input$n.omp.threads2)
    }
  })


  output$plot1 <- renderPlot({
    if(input$select == "Generate data"){
      plot(apply(p.s.1()$p.w.0, 1, mean), w.ho(), main = "Predictions for Holdout Data",
           ylab = "w.ho", xlab = "random effect posterior predictive sample locations")
    }else if(input$select == "Upload dataset") {
      plot(apply(p.s.2()$p.w.0, 1, mean), w.ho.upload(), main = "Predictions for Holdout Data",
           ylab = "w.ho", xlab = "random effect posterior predictive sample locations")
    }
  })







  m.r.1 <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spNNGP(y.new.1()~x.new.1()-1, coords=coords.new.1(), starting=starting, method='response', n.neighbors=input$n.neighbors,
             tuning=tuning, priors=priors, cov.model=input$cov.model,
             n.samples=input$n.samples, n.omp.threads=input$n.omp.threads, n.report=input$n.report)
    }
  })

  m.r.2 <- eventReactive(input$go2, {
    req(data())
    req(input$select)
    if(input$select == "Upload dataset"){
      spNNGP(y.new.2()~x.new.2()-1, coords=coords.new.2(), starting=starting, method='response', n.neighbors=input$n.neighbors2,
             tuning=tuning, priors=priors, cov.model=input$cov.model2,
             n.samples=input$n.samples2, n.omp.threads=input$n.omp.threads2, n.report=input$n.report2)
    }
  })

  p.r.1 <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spNNGP::spPredict(m.r.1(), X.0 = x.ho(), coords.0 = coords.ho(), n.omp.threads=input$n.omp.threads)
    }
  })

  p.r.2 <- eventReactive(input$go2, {
    req(data())
    req(input$select)
    if(input$select == "Upload dataset"){
      spNNGP::spPredict(m.r.2(), X.0 = x.ho.upload(), coords.0 = coords.ho.upload(), n.omp.threads=input$n.omp.threads2)
    }
  })

  output$plot2 <- renderPlot({
    if(input$select == "Generate data"){
      plot(apply(p.s.1()$p.y.0, 1, mean), y.ho(), main = "Predictions for Holdout Data",
           ylab = "y.ho", xlab = "response variable posterior predictive sample locations", col = "red", pch = 17)
      points(apply(p.r.1()$p.y.0, 1, mean), y.ho(), pch=19, col="blue")
      legend("bottomright", legend = c("Sequential Model", "Response Model"), col = c("red", "blue"), pch = c(17, 19))
    }else if(input$select == "Upload dataset"){
      plot(apply(p.s.2()$p.y.0, 1, mean), y.ho.upload(), main = "Predictions for Holdout Data",
           ylab = "y.ho", xlab = "response variable posterior predictive sample locations", col = "red", pch = 17)
      points(apply(p.r.2()$p.y.0, 1, mean), y.ho.upload(), pch=19, col="blue")
      legend("bottomright", legend = c("Sequential Model", "Response Model"), col = c("red", "blue"), pch = c(17, 19))
    }

  })

}


shinyApp(ui = ui, server = server)





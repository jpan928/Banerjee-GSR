rm(list = ls())

library(shiny)
library(spBayes)
library(DT)

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                 "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
verbose <- TRUE


ui <- fluidPage(
  pageWithSidebar(
    headerPanel('spBayes spLM Function: Fittings of Gaussian Univariate Bayesian Spatial Regression Models'),

    sidebarPanel(
      selectInput(inputId = "select", label = "Choose whether you wish to generate data or upload your own dataset",
                  choices = c("Generate data", "Upload dataset"), selected = NULL),

    tags$hr(),

    conditionalPanel(

      condition = "input.select == 'Generate data'",

      numericInput(inputId = "n", 'Number of observations per sample', value = 100),
      numericInput(inputId = "nsamp", "Number of Samples", value = 200),
      numericInput(inputId = "sig.sq", "Sigma Squared", value = 2),
      numericInput(inputId = "tau.sq", "Tau Squared", value = 0.1),
      selectInput(inputId = 'cov.model', "Covariance Model", choices = c('exponential', 'matern', 'spherical', 'gaussian'), selected = 'exponential',
                  multiple = FALSE),
      textInput(inputId = "B", "B matrix: Enter a vector (comma delimited)", value = "1, 5"),
      numericInput(inputId = "phi", "Phi", value = 3 / 0.5),
      numericInput(inputId = 'n.report', 'Interval to report sampling progress', value = 500),
      #radioButtons("verbose", "Verbose (T/F)", choices = c(TRUE, FALSE), selected = TRUE),
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

      numericInput(inputId = "n2", 'Number of observations per sample', value = 100),
      numericInput(inputId = "nsamp2", "Number of Samples", value = 200),
      numericInput(inputId = "sig.sq2", "Sigma Squared", value = 2),
      numericInput(inputId = "tau.sq2", "Tau Squared", value = 0.1),
      selectInput(inputId = 'cov.model2', "Covariance Model", choices = c('exponential', 'matern', 'spherical', 'gaussian'), selected = 'exponential',
                  multiple = FALSE),
      textInput(inputId = "B2", "B matrix: Enter a vector (comma delimited)", value = "1, 5"),
      numericInput(inputId = "phi2", "Phi", value = 3 / 0.5),
      numericInput(inputId = 'n.report2', 'Interval to report sampling progress', value = 500),
      #radioButtons("verbose2", "Verbose (T/F)", choices = c(TRUE, FALSE), selected = TRUE),
      actionButton("go2", "Submit")
    )
    ),

    mainPanel(
      tabsetPanel( type = "tabs",
                   tabPanel("Plot", plotOutput('Plot1')),
                   tabPanel("Quantiles of Theta and Beta Recover Samples (priors.1)", textOutput('label1'),
                            dataTableOutput('txt1'),
                            dataTableOutput('txt2')),
                   tabPanel("Quantiles of Theta and Beta Recover Samples (priors.2)", textOutput('label2'),
                            dataTableOutput('txt3'),
                            dataTableOutput('txt4')),
                   tabPanel("Quantiles of Theta and Beta Predictive Samples (priors.1)", textOutput('label3'),
                            dataTableOutput('txt5'),
                            dataTableOutput('txt6')),
                   tabPanel("Quantiles of Theta and Beta Predictive Samples (priors.2)", textOutput('label4'),
                            dataTableOutput('txt7'),
                            dataTableOutput('txt8'))

                  )

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

  ####################
  # Generate Dataset #
  ####################

  p <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      length(as.numeric(unlist(strsplit(input$B, ","))))
    }
  })


  priors.1.gen <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      list("beta.Norm"=list(rep(0,p()), diag(1000,p())),
           "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
           "tau.sq.IG"=c(2, 0.1))
    }
  })

  coords <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      cbind(runif(input$n,0,1), runif(input$n,0,1))
    }
  })

  X <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      as.matrix(cbind(1, rnorm(input$n)))
    }
  })

  D <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      as.matrix(dist(coords()))
    }
  })

  R <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      exp(-input$phi*D())
    }
  })

  w <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      rmvn(1, rep(0,input$n), input$sig.sq*R())
    }
  })

  y <- reactive({
    req(input$select)
    if(input$select == "Generate data"){
      rnorm(input$n, X()%*%as.numeric(unlist(strsplit(input$B, ","))) + w(), sqrt(input$tau.sq))
    }
  })



  ##################
  # Upload Dataset #
  ##################
  data <- reactive({
    req(input$file1)
    req(input$select)
    if(input$select == "Upload dataset"){
      read.csv(input$file1$datapath, header = input$header, sep = input$sep)
    }
  })

  p2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      length(as.numeric(unlist(strsplit(input$B2, ","))))
    }
  })

  priors.1.up <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      list("beta.Norm"=list(rep(0,p2()), diag(1000,p2())),
           "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
           "tau.sq.IG"=c(2, 0.1))
    }
  })

  coords2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      cbind(data()[,which(colnames(data()) == "coords1")],
            data()[,which(colnames(data()) == "coords2")])
    }
  })

  X2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      cbind(1, data()[,which(colnames(data()) == "x")])
    }
  })

  D2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      as.matrix(dist(coords2()))
    }
  })

  R2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      exp(-input$phi2*D2())
    }
  })

  w2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      rmvn(1, rep(0, nrow(data())), input$sig.sq2*R2())
    }
  })

  y2 <- eventReactive(data(), {
    req(input$select)
    if(input$select == "Upload dataset"){
      rnorm(nrow(data()), X2()%*%as.numeric(unlist(strsplit(input$B2, ","))) + w2(), sqrt(input$tau.sq2))
    }
  })



  m.1.gen <- eventReactive(input$go,{
    req(input$select)
    if(input$select == "Generate data"){
    spBayes::spLM(y()~X()-1, coords=coords(), starting=starting,
                        tuning=tuning, priors=priors.1.gen(), cov.model=input$cov.model,
                        n.samples=input$nsamp, verbose=verbose, n.report=input$n.report)
    }
    })
  m.1.up <- eventReactive(input$go2,{
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spBayes::spLM(y2()~X2()-1, coords=coords2(), starting=starting,
           tuning=tuning, priors=priors.1.up(), cov.model=input$cov.model2,
           n.samples=input$nsamp2, verbose=verbose, n.report=input$n.report2)
    }
  })



  m.2.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spLM(y()~X()-1, coords=coords(), starting=starting,
           tuning=tuning, priors=priors.2, cov.model=input$cov.model,
           n.samples=input$nsamp, verbose=verbose, n.report=input$n.report)
    }
  })
  m.2.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spLM(y2()~X2()-1, coords=coords2(), starting=starting,
           tuning=tuning, priors=priors.2, cov.model=input$cov.model2,
           n.samples=input$nsamp2, verbose=verbose, n.report=input$n.report2)
    }
  })


  burn.in.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      0.5*input$nsamp
    }
  })
  burn.in.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      0.5*input$nsamp2
    }
  })



  ###########################################
  ##recover beta and spatial random effects##
  ###########################################
  m.1.new.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spRecover(m.1.gen(), start=burn.in.gen(), verbose=FALSE)
    }
  })
  m.1.new.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spRecover(m.1.up(), start=burn.in.up(), verbose=FALSE)
    }
  })


  m.2.new.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spRecover(m.2.gen(), start=burn.in.gen(), verbose=FALSE)
    }
  })
  m.2.new.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spRecover(m.2.up(), start=burn.in.up(), verbose=FALSE)
    }
  })


  output$label1 <- renderText({"The first table displays the quantile values of the recovered beta samples
    using the first set of priors (priors.1). The second table displays the quantile values of
    the recovered theta samples also using the first set of priors."})
  output$txt1 <- renderDataTable({
    if(input$select == "Generate data"){
      round(summary(m.1.new.gen()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
      #round(summary(m.2.new.gen()$p.theta.recover.samples)$quantiles[,c(3,1,5],2)
    }else if(input$select == "Upload dataset"){
      round(summary(m.1.new.up()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
      #round(summary(m.2.new.up()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))

  output$txt2 <- renderDataTable({
    if(input$select == "Generate data"){
      round(summary(m.1.new.gen()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
      #round(summary(m.2.new.gen()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      round(summary(m.1.new.up()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
      #round(summary(m.2.new.up()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))



  output$label2 <- renderText({"The first table displays the quantile values of the recovered beta samples
    using the second set of priors (priors.2). The second table displays the quantile values of
    the recovered theta samples also using the second set of priors."})
  output$txt3 <- renderDataTable({
    if(input$select == "Generate data"){
      #round(summary(m.1.new.gen()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
      round(summary(m.2.new.gen()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      #round(summary(m.1.new.up()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
      round(summary(m.2.new.up()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))

  output$txt4 <- renderDataTable({
    if(input$select == "Generate data"){
      #round(summary(m.1.new.gen()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
      round(summary(m.2.new.gen()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      #round(summary(m.1.new.up()$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
      round(summary(m.2.new.up()$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))







  m.1.w.summary.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      summary(mcmc(t(m.1.new.gen()$p.w.recover.samples)))$quantiles[,c(3,1,5)]
    }
  })
  m.1.w.summary.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      summary(mcmc(t(m.1.new.up()$p.w.recover.samples)))$quantiles[,c(3,1,5)]
    }
  })

  m.2.w.summary.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      summary(mcmc(t(m.2.new.gen()$p.w.recover.samples)))$quantiles[,c(3,1,5)]
    }
  })
  m.2.w.summary.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      summary(mcmc(t(m.2.new.up()$p.w.recover.samples)))$quantiles[,c(3,1,5)]
    }
  })



  output$Plot1 <- renderPlot({
    if(input$select == "Generate data"){
      plot(w(), m.1.w.summary.gen()[,1], xlab="Observed w", ylab="Fitted w",
           xlim=range(w()), ylim=range(m.1.w.summary.gen()), main="Spatial random effects")
      arrows(w(), m.1.w.summary.gen()[,1], w(), m.1.w.summary.gen()[,2], length=0.02, angle=90)
      arrows(w(), m.1.w.summary.gen()[,1], w(), m.1.w.summary.gen()[,3], length=0.02, angle=90)
      lines(range(w()), range(w()))

      points(w(), m.2.w.summary.gen()[,1], col="blue", pch=19, cex=0.5)
      arrows(w(), m.2.w.summary.gen()[,1], w(), col="blue", m.2.w.summary.gen()[,2], length=0.02, angle=90)
      arrows(w(), m.2.w.summary.gen()[,1], w(), col="blue", m.2.w.summary.gen()[,3], length=0.02, angle=90)
    }else if(input$select == "Upload dataset"){
      plot(w2(), m.1.w.summary.up()[,1], xlab="Observed w", ylab="Fitted w",
           xlim=range(w2()), ylim=range(m.1.w.summary.up()), main="Spatial random effects")
      arrows(w2(), m.1.w.summary.up()[,1], w2(), m.1.w.summary.up()[,2], length=0.02, angle=90)
      arrows(w2(), m.1.w.summary.up()[,1], w2(), m.1.w.summary.up()[,3], length=0.02, angle=90)
      lines(range(w2()), range(w2()))

      points(w2(), m.2.w.summary.up()[,1], col="blue", pch=19, cex=0.5)
      arrows(w2(), m.2.w.summary.up()[,1], w2(), col="blue", m.2.w.summary.up()[,2], length=0.02, angle=90)
      arrows(w2(), m.2.w.summary.up()[,1], w2(), col="blue", m.2.w.summary.up()[,3], length=0.02, angle=90)
    }
    })




  ############################
  ##Predictive process model##
  ############################

  m.1.pred.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spLM(y()~X()-1, coords=coords(), knots=c(6,6,0.1), starting=starting,
           tuning=tuning, priors=priors.1.gen(), cov.model=input$cov.model,
           n.samples=input$nsamp, verbose=verbose, n.report=input$n.report)
    }
  })
  m.1.pred.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spLM(y2()~X2()-1, coords=coords2(), knots=c(6,6,0.1), starting=starting,
           tuning=tuning, priors=priors.1.up(), cov.model=input$cov.model2,
           n.samples=input$nsamp2, verbose=verbose, n.report=input$n.report2)
    }
  })


  m.2.pred.gen <- eventReactive(input$go, {
    req(input$select)
    if(input$select == "Generate data"){
      spLM(y()~X()-1, coords=coords(), knots=c(6,6,0.1), starting=starting,
           tuning=tuning, priors=priors.2, cov.model=input$cov.model,
           n.samples=input$nsamp, verbose=verbose, n.report=input$n.report)
    }
  })
  m.2.pred.up <- eventReactive(input$go2, {
    req(input$select)
    req(data())
    if(input$select == "Upload dataset"){
      spLM(y2()~X2()-1, coords=coords2(), knots=c(6,6,0.1), starting=starting,
           tuning=tuning, priors=priors.2, cov.model=input$cov.model2,
           n.samples=input$nsamp2, verbose=verbose, n.report=input$n.report2)
    }
  })



  output$label3 <- renderText({"The first table displays the quantile values of the predicted beta samples
    using the first set of priors (priors.1). The second table displays the quantile values of
    the predicted theta samples also using the first set of priors."})
  output$txt5 <- renderDataTable({
    if(input$select == "Generate data"){
      round(summary(window(m.1.pred.gen()$p.beta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
      #round(summary(window(m.2.pred.gen()$p.beta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      round(summary(window(m.1.pred.up()$p.beta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
      #round(summary(window(m.2.pred.up()$p.beta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))

  output$txt6 <- renderDataTable({
    if(input$select == "Generate data"){
      round(summary(window(m.1.pred.gen()$p.theta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
      #round(summary(window(m.2.pred.gen()$p.beta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      round(summary(window(m.1.pred.up()$p.theta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
      #round(summary(window(m.2.pred.up()$p.beta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))


  output$label4 <- renderText({"The first table displays the quantile values of the predicted beta samples
    using the second set of priors (priors.2). The second table displays the quantile values of
    the predicted theta samples also using the second set of priors."})

  output$txt7 <- renderDataTable({
    if(input$select == "Generate data"){
      #round(summary(window(m.1.pred.gen()$p.theta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
      round(summary(window(m.2.pred.gen()$p.beta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      #round(summary(window(m.1.pred.up()$p.theta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
      round(summary(window(m.2.pred.up()$p.beta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))

  output$txt8 <- renderDataTable({
    if(input$select == "Generate data"){
      #round(summary(window(m.1.pred.gen()$p.theta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
      round(summary(window(m.2.pred.gen()$p.theta.samples, start=burn.in.gen()))$quantiles[,c(3,1,5)],2)
    }else if(input$select == "Upload dataset"){
      #round(summary(window(m.1.pred.up()$p.theta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
      round(summary(window(m.2.pred.up()$p.theta.samples, start=burn.in.up()))$quantiles[,c(3,1,5)],2)
    }
  }, options = list(include.rownames = TRUE, initComplete = JS(
    "function(settings, json) {",
    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
    "}")))

}

shinyApp(ui = ui, server = server)


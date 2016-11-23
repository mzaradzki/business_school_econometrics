library(shiny)
library(MASS)

# Define server logic for random distribution application
shinyServer(function(input, output) {
  
  getIntercept <- reactive({  
    as.numeric(input$intercept)
  })
  
  getCoefficient <- reactive({  
    as.numeric(input$coefficient)
  })
  
  getStdDeviation <- reactive({  
    as.numeric(input$stdeviation)
  })
  
  getNbObs <- reactive({  
    as.numeric(input$nbobs)
  })
  
  # number of explanatory variables
  K <- reactive({  
    1
  })
  
  # degree of freedom
  DoF <- reactive({  
    getNbObs() - K() - 1
  })
  
  Xs <- reactive({  
    seq(1, getNbObs())
  })
  
  Es <- reactive({  
    errors <- rnorm(getNbObs(), 0, getStdDeviation())
    (errors - mean(errors))
    #errors
  })
  
  Ys <- reactive({  
    getIntercept() + getCoefficient()*Xs() + Es()
  })
  
  DF <- reactive({  
    data.frame(X=Xs(), X2=Xs()*Xs(), X3=Xs()*Xs()*Xs(), E=Es(), Y=Ys())
  })
  
  FIT <- reactive({  
    lm(Y ~ X, DF())
  })
  
  FIT123 <- reactive({  
    lm(Y ~ X + X2 + X3, DF())
  })
  
  R2 <- reactive({  
    summary(FIT())$r.squared
  })
  
  EMPFISH <- reactive({  
    empfish <- (R2()/(1-R2())) * (DoF()/K())
  })
  
  output$summary <- renderPrint({
    summary(FIT(), correlation = TRUE)[c('call','coefficients','sigma','r.squared','correlation')]
  })
  
  FITS <- reactive({  
    M = matrix(0,100,2)
    for (i in 1:100) {
      locE <- rnorm(getNbObs(), 0, getStdDeviation())
      locY <- getIntercept() + getCoefficient()*Xs() + locE
      locDF <- data.frame(X=Xs(), Y=locY)
      locFIT <- lm(Y ~ X, locDF)
      locBETAS = matrix(locFIT$coefficients)
      M[i,1] <- locBETAS[1]
      M[i,2] <- locBETAS[2]
    }
    M
  })
  
  output$plotFit <- renderPlot({  
    LOW <- min(min(Ys()), min(FIT()$fitted.values))
    HIGH <- max(max(Ys()), max(FIT()$fitted.values))
    YRANGE <- c(LOW, HIGH)
    plot(x=Ys(), y=FIT()$fitted.values, ylab="Fitted Value", xlab="Original Value", col="red", type="p", ylim=YRANGE)
    lines(x=Ys(), y=Ys(), col="black", type="l", ylim=YRANGE)
  })
  
  output$plotError <- renderPlot({  
    hist(FIT()$residual)
  })
  
  output$plotBoxes <- renderPlot({  
    par(mfrow=c(1,2))
    boxplot(FITS()[,1], ylab="estimates", main="Intercept")
    boxplot(FITS()[,2], ylab="estimates", main="Beta")
  })
  
  output$plotBoxIntercept <- renderPlot({  
    boxplot(FITS()[,1], ylab="intercept")
  })
  
  output$plotBoxBeta <- renderPlot({  
    boxplot(FITS()[,2], ylab="coefficient")
  })
  
  studentChart <- function(tratio){
    x <- seq(-6, 6, length=500)
    y <- dt(x, DoF(), log=FALSE)
    plot(x, y, ylab="probability density", type="l", lwd=1)
    if(tratio > 0) {
      if(max(x) >= tratio) {
        ind1 <- min(which(x >= tratio))
        ind2 <- length(x)
        polygon(x=c(x[c(ind1, ind1:ind2, ind2)]), y=c(0, y[ind1:ind2], 0),col='skyblue')
      }
    }
    if(tratio < 0) {
      if(min(x) <= tratio) {
        ind1 <- 1
        ind2 <- max(which(x <= tratio))
        polygon(x=c(x[c(ind1, ind1:ind2, ind2)]), y=c(0, y[ind1:ind2], 0),col='skyblue')
      }
    }
  }
  
  output$plotStudentIntercept <- renderPlot({  
    betameans <- FIT()$coefficients
    betastdevs <- sqrt(diag(vcov(FIT())))
    studentChart(betameans[1]/betastdevs[1])
  })
  
  output$plotStudentBeta <- renderPlot({  
    betameans <- FIT()$coefficients
    betastdevs <- sqrt(diag(vcov(FIT())))
    studentChart(betameans[2]/betastdevs[2])
  })
  
  fisherChart <- function(dof1, dof2, empfish){
    top <- max(empfish*1.2, qf(0.99, dof1, dof2, log=FALSE))
    x <- seq(0, top, length=500)
    y <- pf(x, dof1, dof2, log=FALSE)
    plot(x, y, ylab="cummulative probability", type="l", lwd=1)
    if(max(x) >= empfish) {
      ind1 <- max(which(x <= empfish))
      ind2 <- min(which(x >= empfish))
      polygon(x=c(x[c(ind1, ind1:ind2, ind2)]), y=c(0, y[ind1:ind2], 0),col='red')
    }
  }
  
  output$plotFisherR2 <- renderPlot({  
    fisherChart(K(), DoF(), EMPFISH())
  })
  
  output$plotBiNorm <- renderPlot({  
    ##set.seed(17)
    
    center <- FIT()$coefficients
    sigma <- vcov(FIT())
    sigma.inv = solve(sigma)
    ellipse <- function(s,t) {u<-c(s,t)-center; u %*% sigma.inv %*% u / 2}
    
    p <- mvrnorm(1000, center, sigma)
    
    n <- 50
    x <- seq(min(p[,1]), max(p[,1]), length=n)
    y <- seq(min(p[,2]), max(p[,2]), length=n)
    z <- mapply(ellipse, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, "+")))
    plot(t(center), pch=20, xlim=c(min(p[,1]),max(p[,1])), ylim=c(min(p[,2]),max(p[,2])), xlab="Intercept", ylab="Beta")
    contour(x, y, matrix(z,n,n), levels=(0:10), col = terrain.colors(11), add=TRUE)
  })
  
  output$plotFisherJoint123 <- renderPlot({  
    # for a square+cube model
    p2 <- K()+2+1
    p1 <- K()+1
    rss2 <- sum(FIT123()$residual * FIT123()$residual)
    rss1 <- sum(FIT()$residual * FIT()$residual)
    #
    jntfsh <- ((rss1 - rss2)/(p2-p1)) / (rss2/(getNbObs()-p2))
    fisherChart(p2-p1, getNbObs()-p2, jntfsh)
  })
  
  output$plotFisherJointTrue <- renderPlot({  
    # for the real model
    p2 <- K()+1
    p1 <- 0
    rss2 <- sum(FIT()$residual * FIT()$residual)
    rss1 <- sum(Es() * Es())
    #
    jntfsh <- ((rss1 - rss2)/(p2-p1)) / (rss2/(getNbObs()-p2))
    fisherChart(p2-p1, getNbObs()-p2, jntfsh)
  })
  
  chisquareChart <- function(){
    low <- qchisq(0.001, DoF())
    top <- qchisq(0.999, DoF())
    x <- seq(low, top, length=500)
    y <- dchisq(x, df = DoF())
    plot((x/DoF())-1, y*DoF(), xlab="(VarEst-Var)/Var", ylab="probability density", type="l", lwd=1)
    # lowest 5 percentile
    indQLow1 <- 1
    indQLow2 <- max(which(x <= qchisq(0.05, DoF())))
    polygon(x=c(x[c(indQLow1,indQLow1:indQLow2,indQLow2)]/DoF())-1, y=c(0,y[indQLow1:indQLow2],0)*DoF(),col='skyblue')
    #
    # highest 5 percentile
    indQTop1 <- min(which(x >= qchisq(0.95, DoF())))
    indQTop2 <- length(x)
    polygon(x=c(x[c(indQTop1,indQTop1:indQTop2,indQTop2)]/DoF())-1, y=c(0,y[indQTop1:indQTop2],0)*DoF(),col='skyblue')
  }
  
  output$plotRelativeSigmaError <- renderPlot({  
    chisquareChart()
  })
  
})
library(shiny)
library(Hmisc)

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
  
  Xs <- reactive({  
    #seq(1, getNbObs())
    #rnorm(getNbObs(), 0, 1)
    runif(getNbObs(),-10,10)
  })
  
  Es <- reactive({  
    rnorm(getNbObs(), 0, getStdDeviation())
  })
  
  Ys <- reactive({  
    getIntercept() + getCoefficient()*Xs() + Es()
  })
  
  DF <- reactive({  
    data.frame(X=Xs(), E=Es(), Y=Ys())
  })
  
  FIT <- reactive({  
    lm(Y ~ X, DF())
  })
  
  getNbSims <- reactive({  
    100
  })
  
  FITS <- reactive({  
    M = matrix(0,getNbSims(),4)
    for (i in 1:getNbSims()) {
      locE <- rnorm(getNbObs(), 0, getStdDeviation())
      locY <- getIntercept() + getCoefficient()*Xs() + locE
      locDF <- data.frame(X=Xs(), Y=locY)
      locFIT <- lm(Y ~ X, locDF)
      locBETAS = matrix(locFIT$coefficients)
      M[i,1] <- locBETAS[1]
      M[i,2] <- locBETAS[2]
      M[i,3] <- summary(locFIT)$sigma
      M[i,4] <- summary(locFIT)$r.squared
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
    par(mfrow=c(2,2))
    boxplot(FITS()[,1], ylab="estimates", main="Intercept")
    boxplot(FITS()[,2], ylab="estimates", main="Beta")
    boxplot(FITS()[,3], ylab="estimates", main="Sigma")
    boxplot(FITS()[,4], ylab="estimates", main="R2")
  })
  
  output$plotBoxIntercept <- renderPlot({  
    boxplot(FITS()[,1], ylab="intercept")
  })
  
  output$plotBoxBeta <- renderPlot({  
    boxplot(FITS()[,2], ylab="coefficient")
  })
  
  normChart <- function(avg, stdev, title){
    x <- seq(avg-5*stdev, avg+5*stdev, length=500)
    y <- dnorm(x, mean=avg, sd=stdev)
    plot(x, y, type="l", lwd=1, main=title)
    if(avg < 0) {
      if(max(x) >= 0) {
        ind1 <- min(which(x >= 0))
        ind2 <- length(x)
        polygon(x=c(x[c(ind1, ind1:ind2, ind2)]), y=c(0, y[ind1:ind2], 0), col='skyblue')
      }
    }
    if(avg > 0) {
      if(min(x) <= 0) {
        ind1 <- 1
        ind2 <- max(which(x <= 0))
        polygon(x=c(x[c(ind1, ind1:ind2, ind2)]), y=c(0, y[ind1:ind2], 0), col='skyblue')
      }
    }
  }
  
  output$plotNormIntercept <- renderPlot({  
    betameans <- FIT()$coefficients
    betastdevs <- sqrt(diag(vcov(FIT())))
    normChart(betameans[1], betastdevs[1], "Intercept")
  })
  
  output$plotNormBeta <- renderPlot({  
    betameans <- FIT()$coefficients
    betastdevs <- sqrt(diag(vcov(FIT())))
    normChart(betameans[2], betastdevs[2], "Beta")
  })
  
  output$plotNorms <- renderPlot({  
    par(mfrow=c(1,2))
    betameans <- FIT()$coefficients
    betastdevs <- sqrt(diag(vcov(FIT())))
    normChart(betameans[1], betastdevs[1], "Intercept")
    normChart(betameans[2], betastdevs[2], "Beta")
  })
  
  # Generate a summary of the data
  output$summary <- renderPrint({
    summary(FIT())[c('call','coefficients','sigma','r.squared')]
  })
  
  getNewX <- reactive({  
    as.numeric(input$newX)
  })
  
  output$plotHistForecast <- renderPlot({  
    lows <- FITS()[,1] + FITS()[,2]*getNewX() - 4*FITS()[,3]
    tops <- FITS()[,1] + FITS()[,2]*getNewX() + 4*FITS()[,3]
    #
    Yrange = seq(min(lows), max(tops), length=100)
    D = seq(0, 0, length=100)
    #
    for (i in 1:getNbSims()) {
      component <- dnorm(Yrange, mean=(FITS()[i,1]+FITS()[i,2]*getNewX()), sd=FITS()[i,3], log=FALSE)
      D <- D+component
    }
    #
    Ymean <- wtd.mean(Yrange, D, normwt=TRUE)
    Yvar <- wtd.var(Yrange, D, normwt=TRUE)
    #
    title <- paste("Forecast with X=",as.character(getNewX())," (mean=",format(Ymean,digits=1),", stdev=",format(sqrt(Yvar),digits=1),")",sep="")
    #
    plot(x=Yrange, y=D, xlab="Y values", ylab="probability density", main=title, type="l", col="black")
  })
  
})
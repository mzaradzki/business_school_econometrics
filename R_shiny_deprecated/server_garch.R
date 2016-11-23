library(shiny)
library(MASS)
library(Quandl)
library(tseries)
library(moments)
library(fGarch)
library(lmtest)




Quandl.auth("8aVZ3CksqJmfEozdFNbB")

FF = Quandl("KFRENCH/FACTORS_W", trim_start="2000-01-01", trim_end="2014-06-15")
colnames(FF) <- c("Date", "NETMKT", "SMB", "HML", "RF")
FF <- transform(FF, Date=rev(Date), NETMKT=rev(NETMKT), SMB=rev(SMB), HML=rev(HML))

# The yield associated to the RF weekly P&L
FF$RFY <- FF$RF*52/100

SPDR = Quandl(c("GOOG/AMEX_XLF.4", "GOOG/AMEX_XLY.4", "GOOG/AMEX_XLV.4", "GOOG/AMEX_XLI.4", "GOOG/AMEX_XLK.4", "GOOG/AMEX_XLB.4", "YAHOO/INDEX_VIX.CLOSE"), trim_start="2000-01-01", trim_end="2014-06-15")
colnames(SPDR) <- c("Date", "XLF", "XLY", "XLV", "XLI", "XLK", "XLB", "VIX")

SPDR <- transform(SPDR, Date=rev(Date), XLF=rev(XLF), XLY=rev(XLY), XLV=rev(XLV), XLI=rev(XLI), XLK=rev(XLK), XLB=rev(XLB), VIX=rev(VIX))

CDATES <- intersect(SPDR$Date, FF$Date)

##which(SPDR$Date %in% CDATES)

SPDR <- SPDR[which(SPDR$Date %in% CDATES),]

SPDR <- transform(SPDR, XLF=100*log(XLF), XLY=100*log(XLY), XLV=100*log(XLV), XLI=100*log(XLI), XLK=100*log(XLK), XLB=100*log(XLB))

SPDR <- transform(SPDR, DXLF=c(0,diff(XLF,lag=1)), DXLY=c(0,diff(XLY,lag=1)), DXLV=c(0,diff(XLV,lag=1)), DXLI=c(0,diff(XLI,lag=1)), DXLK=c(0,diff(XLK,lag=1)), DXLB=c(0,diff(XLB,lag=1)))

FF <- FF[which(FF$Date %in% CDATES),]

ALL <- merge(FF, SPDR, by="Date")

# Define server logic for random distribution application
shinyServer(function(input, output) {
  
  
  
  DF <- reactive({  
    ## can we make a copy ?
    ALL
  })
  
  getSector <- reactive({  
    switch(input$sector,
              XLV = 'XLV',
              XLF = 'XLF',
              XLI = 'XLI',
              XLK = 'XLK',
              XLY = 'XLY',
              XLB = 'XLB',
              XLV)
  })
  
  getCoefficient <- reactive({  
    as.numeric(input$coefficient)
  })
  
  getStdDeviation <- reactive({  
    as.numeric(input$stdeviation)
  })
  
  getNbObs <- reactive({  
    length(DF()$Date)
  })
  
  # number of explanatory variables
  K <- reactive({  
    3
  })
  
  # degree of freedom
  DoF <- reactive({  
    getNbObs() - K() - 1
  })
  
  FITfama <- reactive({  
    localDF <- DF()
    localDF$DSECTOR <- matrix(t(localDF[paste("D", getSector(), sep="")]))
    lm(DSECTOR ~ NETMKT + SMB + HML, localDF, y=TRUE)
  })
  
  FITcapm <- reactive({  
    localDF <- DF()
    localDF$DSECTOR <- matrix(t(localDF[paste("D", getSector(), sep="")]))
    lm(DSECTOR ~ NETMKT, localDF, y=TRUE)
  })
  
  ERRORS <- reactive({  
    XLFerr <- lm(DXLF ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    XLVerr <- lm(DXLV ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    XLYerr <- lm(DXLY ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    XLKerr <- lm(DXLK ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    XLIerr <- lm(DXLI ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    XLBerr <- lm(DXLB ~ NETMKT + SMB + HML, DF(), y=TRUE)$residual
    data.frame(XLF=XLFerr, XLV=XLVerr, XLY=XLYerr, XLK=XLKerr, XLI=XLIerr, XLB=XLBerr)
  })
  
  DW <- reactive({  
    dwtest(FITfama())
  })
  
  BG <- reactive({  
    bgtest(FITfama(), order=2, type="F")
  })
  
  R2 <- reactive({  
    summary(FITfama())$r.squared
  })
  
  EMPFISH <- reactive({  
    empfish <- (R2()/(1-R2())) * (DoF()/K())
  })
  
  output$summary <- renderPrint({
    summary(FITfama(), correlation = TRUE)[c('call','coefficients','sigma','r.squared','correlation')]
  })
  
  output$durbinwatson <- renderPrint({
    DW()
  })
  
  output$breuschgodfrey <- renderPrint({
    BG()
  })
  
  GARCH <- reactive({  
    p <- as.integer(as.numeric(input$garchP))
    q <- as.integer(as.numeric(input$garchQ))
    ##fit <- garchFit(~garch(p, q), data = DF()$XLF)
    fit <- garchFit(~garch(1, 1), data = FITfama()$residual, trace=FALSE)
    fit
  })
  
  output$agostinoskewness <- renderPrint({
    agostino.test(FITfama()$residual)[c('statistic','p.value')]
  })
  
  output$jarqueberranormality <- renderPrint({
    jarque.bera.test(FITfama()$residual)
  })
  
  output$plotFit <- renderPlot({  
    LOW <- min(min(FITfama()$y), min(FITfama()$fitted.values))
    HIGH <- max(max(FITfama()$y), max(FITfama()$fitted.values))
    YRANGE <- c(LOW, HIGH)
    plot(x=FITfama()$y, y=FITfama()$fitted.values, ylab="Fitted Value", xlab="Original Value", col="red", type="p", ylim=YRANGE)
    lines(x=FITfama()$y, y=FITfama()$y, col="black", type="l", ylim=YRANGE)
  })
  
  output$plotError <- renderPlot({  
    par(mfrow=c(1,2))
    plot(FITfama()$residual, main="Time-series", type="l")
    ##
    dindex <- FITfama()$residual
    ##dindex.garch <- garch(dindex, order = c(1,1), trace=FALSE)
    y1 <- dindex
    y1 <- y1[!is.na(y1)]
    ##y1 <- y1 - mean(y1)
    ##y1 <- y1 / sqrt(var(y1))
    ##y2 <- dindex.garch$residuals
    ##y2 <- y2[!is.na(y2)]
    ##y2 <- y2 - mean(y2)
    ##y2 <- y2 / sqrt(var(y2))
    qqnorm(y1, col="red")
    ##par( new = TRUE )
    ##qqnorm(y2, xlim=c(-4, +4), ylim=c(-4, +4), col="blue")
  })
  
  output$plotScaledError <- renderPlot({  
    plot(FITfama()$residual/DF()$VIX, type="l")
  })
  
  output$plotQQError <- renderPlot({  
    dindex <- FITfama()$residual
    ##dindex.garch <- garch(dindex, order = c(1,1), trace=FALSE)
    y1 <- dindex
    y1 <- y1[!is.na(y1)]
    y1 <- y1 - mean(y1)
    y1 <- y1 / sqrt(var(y1))
    ##y2 <- dindex.garch$residuals
    ##y2 <- y2[!is.na(y2)]
    ##y2 <- y2 - mean(y2)
    ##y2 <- y2 / sqrt(var(y2))
    qqnorm(y1, col="red")
    ##par( new = TRUE )
    ##qqnorm(y2, xlim=c(-4, +4), ylim=c(-4, +4), col="blue")
  })
  
  output$plotGARCHError <- renderPlot({  
    par(mfrow=c(1,2))
    residuals <- FITfama()$residual
    volatilities <- volatility(GARCH(),type="sigma")
    plot(residuals, main="Time-series", col="red", type="l", xlab=NA, ylab=NA)
    par(new = TRUE)
    plot(volatilities, col="black", type="l", axes=FALSE, xlab=NA, ylab=NA)
    ##axis(side = 4)
    ##mtext(side = 4, line = 3, "volatility")
    ##
    y1 <- (residuals/volatilities)
    y1 <- y1[!is.na(y1)]
    y1 <- y1 - mean(y1)
    y1 <- y1 / sqrt(var(y1))
    ##
    qqnorm(y1, col="red")
    ##par( new = TRUE )
    ##qqnorm(y2, xlim=c(-4, +4), ylim=c(-4, +4), col="blue")
  })
  
  
  
  
  
})
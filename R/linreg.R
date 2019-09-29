#' Linreg function for Linear Regression
#' to perform different calculations for Linear Regression
#' @param formula is a formula
#' @param data is a DataFrame name
#' both parameters are inputs for Linear Regression different calculations
#' to find out Regression Coefficients and other factors
#' @return the list of different objects of Linreg class
#' @examples
#' Linreg(Petal.Length~Sepal.Width+Sepal.Length,iris)
#' Linreg(Petal.Length ~ Sepal.Length,iris)
#' @export

linreg <- function(formula,data){
  x <- model.matrix(formula,data)
  y <- as.matrix(data[all.vars(formula)[1]])
  beta_names <- colnames(x)
  call <- match.call()
  #degree of freedom
  n <- nrow(x)
  p <- ncol(x)
  df <- n- p
  beta <- solve(t(x)%*%x) %*% t(x)%*%y
  beta_named <- data.frame(t(beta))
  colnames(beta_named) <- beta_names
  rownames(beta_named) <- " "
  ycap <- c(x%*%beta,dimnames("FittedValues"))
  ecap <- c(y - ycap,dimnames("Residuals"))
  sigm2 <- as.vector((t(ecap) %*% ecap))/df
  varbeta <- sigm2 * as.vector(diag(solve(t(x)%*% x)))
  tvalue <- as.vector(beta)/sqrt(varbeta)
  pvalues <- 2* pt(-abs(tvalue),df)
  a <- list(beta_named,ycap,ecap,df,sigm2,varbeta,tvalue,pvalues,call,beta)
  class(a) <- "linreg"
  names(a) <- c("Regression_Coefficients","Fitted_Values","Residuals",
                "Degree_of_Freedom","Residual_Variance",
                "Variance_of_Regression_Coefficients","tvalues",
                "Pvalues","call","beta_matrix")
  return(a)
}
#Residual Coefficients
coef.linreg <- function(a){
    if(attr(a,"class")=="linreg")
    return(a$Regression_Coefficients)
}
#fitted values
pred.linreg <- function(a){
  if(attr(a,"class")=="linreg")
  return(a$Fitted_Values)
}
#residuals
residuals.linreg <- function(a){
  if(attr(a,"class")=="linreg")
  return(a$Residuals)
}
summary.linreg <- function(a){
  if(attr(a,"class")=="linreg")
    cat("\nCall:\n", paste(deparse(a$call),"\n"),sep ="\n")
    stderrs <- sqrt(a$Variance_of_Regression_Coefficients)
    res_stderr <- sqrt(a$Residual_Variance)
    colnames(a$beta_matrix) <- "Estimate"
    coefficient <- data.frame(Estimate = a$beta_matrix,
                              Std.Error =stderrs,t.value =a$tvalues,p.value = a$Pvalues,estim = "***")
    rownames(coefficient) <- colnames(a$Regression_Coefficients)
    print(coefficient)

    cat("\nResidual standard error:",res_stderr,"on",a$Degree_of_Freedom,"degrees of freedom\n")
}
print.linreg <- function(a){
  if(attr(a,"class")=="linreg"){
    cat("\nCall:\n", paste(deparse(a$call)),sep ="\n")
  }
  if(!is.null(coef(a))){
    cat("\n Coefficients:\n")
    print.data.frame(format(a$Regression_Coefficients),
                  quote = FALSE)
  }
}
plot.linreg <- function(a){
  if(attr(a,"class") == "linreg" ){
    require(ggplot2)
    require(gridExtra)
    require(png)
    require(RCurl)
    #Loading images
    img0 <-readPNG(getURLContent('https://i.imgur.com/ZFdIPwp.png'))
    grob0 <- grid::rasterGrob(img0,height = 1,width = 1,interpolate = TRUE)
    #img <- readPNG("images/LiU_primar_svart.png")
    dfr <- as.data.frame(cbind("Residuals" = a$Residuals,"Fitted_Values" = a$Fitted_Values))
    plot1 <- ggplot(data=dfr,aes(y=Residuals,x=Fitted_Values))+
      geom_point()+
      labs(x="Fitted Values",y="Residuals",title = "Residuals Vs Fitted")+
      theme(plot.background = element_rect(fill = "skyblue",color = "skyblue"))+
      geom_smooth(method = "loess", formula= y~x,se= TRUE)
    dfr2 <- as.data.frame(cbind("Residuals" = a$Residuals,"Fitted_Values" = a$Fitted_Values))
    plot2 <- ggplot(data=dfr,aes(y=sqrt(abs(Residuals)),x=Fitted_Values))+
      geom_point()+
      labs(x="Fitted Values",y=expression(paste(sqrt("Standardized Residuals"))),title = "Scale-Location")+
      theme(plot.background = element_rect(fill = "skyblue",color = "skyblue"))+
      geom_smooth(method = "loess", formula= y~x,se= TRUE)
    grid.arrange(grob0,plot1,plot2)
  }
}
pred <- function(a){
  if(attr(a,"class")=="linreg"){
    prediction <- UseMethod("pred")
    return(prediction)
  }
}

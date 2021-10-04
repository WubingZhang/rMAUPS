#' Peptide/Spectra Count Based Empirical Bayes Statistics for Differential Expression
#'
#' @docType methods
#' @name SpectraCounteBayes
#' @rdname SpectraCounteBayes
#'
#' @param fit an list object produced by Limma \code{eBayes} function, it should
#' have one additional attribute \code{$count}, which stored the peptide or PSM count
#' quantified for the gene in label-free or isobaric labelled data.
#' @param fit.method the method used to fit variance against the number of
#' peptides/PSM count quantified. Two available methods: "loess","nls" and
#' "spline". default "loess"."loess" uses \code{loess} and span = 0.75, "nls""
#' uses a explicit formula \code{y~a+b/x}. "spline" uses \code{smooth.spline} and
#' "generalized cross-validation" for smoothing parameter computation. For "nls",
#' independent variable x is peptide/PSM count, response y is pooled variance
#' (\code{fit$sigma^2}). For "loess" and "spline" method, both x and y are log
#' transformed before applying the two methods. In most of time, "loess" is sufficient.
#' To quickly assess the fit model, use \code{VarianceScatterplot} and
#' \code{Residualplot} functions.
#' @param coef_col an integer vector indicating the column(s) of fit$coefficients
#' for which the function is to be performed. if not specified, all coefficients
#' are used.
#'
#' @return a list object with the additional attributes being: sca.t - Spectra
#' Count Adjusted posterior t-value sca.p - Spectra Count Adjusted posterior
#' p-value sca.dfprior - Spectra Count Adjusted prior degrees of freedom sca.priorvar-
#' Spectra Count Adjusted estimated prior variance sca.postvar - Spectra Count Adjusted
#' posterior variance loess.model - fitted loess model.
#'
#' @details This function adjusts the T-statistics and p-values for quantitative MS
#' proteomics experiment according to the number of peptides/PSMs used for
#' quantification. The method is similar in nature to intensity-based Bayes
#' method (Maureen A. Sartor et al BMC Bioinformatics 2006).
#'
#' @author Yafeng Zhu
#'
#' @export

SpectraCounteBayes<-function(fit,fit.method="loess",coef_col) {

    ################################################
    #  The function was adapted from:
    #  Function for IBMT (Intensity-based Moderated
    #  T-statistic) Written by Maureen Sartor
    #  University of Cincinnati, 2006
    ################################################
    logVAR<-log(fit$sigma^2)
    df<-fit$df.residual
    numgenes<-length(logVAR[df>0])
    df[df==0]<-NA
    eg<-logVAR-digamma(df/2)+log(df/2)
    names(fit$count) <- rownames(fit$coefficients)
    output<-fit
    output$fit.method <- fit.method

    if (fit.method == "loess"){
        x<-log2(fit$count)
        loess.model <- loess(logVAR~x,span = 0.75)
        y.pred <- fitted(loess.model)
        output$model <- loess.model
    }else if (fit.method == "nls"){
        x<-fit$count
        y<-fit$sigma^2
        nls.model <-  nls(y~a+b/x,start = (list(a=0.1,b=0.05)))
        y.pred <- log(fitted(nls.model))
        output$model <- nls.model
    }else if (fit.method == "spline"){
        x<-log2(fit$count)
        spline.model <- smooth.spline(x,logVAR,cv=FALSE)
        y.pred <- fitted(spline.model)
        output$model <- spline.model
    }

    egpred<-y.pred-digamma(df/2)+log(df/2)

    myfct<- (eg-egpred)^2 - trigamma(df/2)

    mean.myfct<-mean(myfct,na.rm=TRUE)

    priordf<-vector()
    testd0<-vector()

    for (i in seq(1,numgenes*10)) {
        testd0[i]<-i/10
        priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
        if (i>2) {
            if (priordf[i-2]<priordf[i-1]) { break }
        }
    }
    d0<-testd0[match(min(priordf),priordf)] # prior degree found

    s02<-exp(egpred + digamma(d0/2) - log(d0/2)) # calculate prior variance

    post.var<- (d0*s02 + df*fit$sigma^2)/(d0+df)
    post.df<-d0+df
    # sca.t and scc.p stands for spectra count adjusted t and p values.
    sca.t<-as.matrix(fit$coefficients[,coef_col]/(fit$stdev.unscaled[,coef_col]
    *sqrt(post.var)))
    sca.p<-as.matrix(2*pt(abs(sca.t),post.df,lower.tail = FALSE))

    output$sca.t<-sca.t
    output$sca.p<-sca.p
    output$sca.postvar<-post.var
    output$sca.priorvar<-s02
    output$sca.dfprior<-d0

    return (output)
}

#' generate a boxplot of the variance
#'
#' @docType methods
#' @name VarianceBoxplot
#' @rdname VarianceBoxplot
#'
#' @param fit an object returned from \code{SpectraCounteBayes} function.
#' @param n set a number to plot only the genes with count value smaller or equal to n.
#' @param xlab the title for x axis.
#' @param ylab the title for y axis.
#' @param main the title for the figure.
#'
#' @return a ggplot instance.
#'
#' @author Yafeng Zhu
#' @import ggplot2 ggpubr
#'
#' @export
VarianceBoxplot <- function (fit, n=20, xlab="PSM count",
                             ylab = "log(Variance)", main=NULL){
  requireNamespace("ggplot2")
  requireNamespace("ggpubr")

  x <- fit$count
  y <- fit$sigma^2

  df.temp <- data.frame("pep_count" =x, "variance" = y )
  df.temp.filter <- df.temp[df.temp$pep_count<=n,]
  df.temp.filter$logVariance = log(df.temp.filter$variance)
  if (fit$fit.method=="nls"){
    y.pred <- log(predict(fit$model,data.frame(x=seq(1,n))))
  }else if (fit$fit.method=="loess"){
    y.pred <- predict(fit$model,data.frame(x=log2(seq(1,n))))}
  else if (fit$fit.method=="spline"){
    y.pred <- predict(fit$model,x=log2(seq(1,n)))$y
  }
  p = ggplot()
  p = p + geom_boxplot(data = df.temp.filter,
                       aes_string("pep_count", "logVariance", group="pep_count"))
  p = p + geom_line(aes_string("pep_count", "y.pred"), col='red',lwd=1,
                    data = na.omit(data.frame(pep_count = seq(1,n),
                                      y.pred = y.pred)))
  p = p + theme_pubr(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5, size=18))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + labs(x = xlab, y = ylab, title = main)
  p
}

#' generate a scatter plot of the variance
#'
#' @docType methods
#' @name VarianceScatterplot
#' @rdname VarianceScatterplot
#'
#' @param fit an object returned from \code{SpectraCounteBayes} function.
#' @param xlab the title for x axis.
#' @param ylab the title for y axis.
#' @param main the title for the figure.
#'
#' @return a ggplot instance.
#'
#' @author Yafeng Zhu
#' @import ggplot2 ggpubr
#'
#' @export
VarianceScatterplot <- function (fit, xlab="log2(count)",
                                 ylab = "log(Variance)",main=NULL){
  x <- fit$count
  y <- fit$sigma^2

  if (fit$fit.method=="nls"){
    y.pred <- log(fitted(fit$model))
  }else if (fit$fit.method=="loess" | fit$fit.method=="spline" ) {
    y.pred <- fitted(fit$model)}

  gg = data.frame(x = log2(x), y = log(y))
  p = ggplot()
  p = p + geom_point(aes(x, y), shape = 1, size = 0.5, data = gg)
  p = p + geom_line(aes(x, y), col='red',lwd=1,
                    data = data.frame(x = log2(x[order(x)]),
                                      y = y.pred[order(x)]))
  p = p + theme_pubr(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5, size=18))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + labs(x = xlab, y = ylab, title = main)
  p
}

#' plot the residuals against the number of quantified peptides/PSMs.
#'
#' @docType methods
#' @name Residualplot
#' @rdname Residualplot
#'
#' @param fit an object returned from \code{SpectraCounteBayes} function.
#' @param xlab the title for x axis.
#' @param ylab the title for y axis.
#' @param main the title for the figure.
#'
#' @return a ggplot instance.
#'
#' @author Yafeng Zhu
#' @import ggplot2 ggpubr
#'
#' @export
Residualplot <- function (fit, xlab="log2(count)",
                          ylab="Variance(fitted - observed)", main=NULL){
  x <- fit$count

  if (fit$fit.method=="nls"){
    y <- log(fitted(fit$model)) - log(fit$sigma^2)
  }else if (fit$fit.method=="loess" | fit$fit.method=="spline") {
    y <- residuals(fit$model)}

  gg = data.frame(x = log2(x), y = y)
  p = ggplot(gg, aes(x, y))
  p = p + geom_point(shape = 20, size = 0.5)
  p = p + theme_pubr(base_size = 14)
  p = p + theme(plot.title = element_text(hjust = 0.5, size=18))
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + labs(x = xlab, y = ylab, title = main)
  p
}


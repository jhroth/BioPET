# Jeremy Roth
library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(VGAM)

enrichment_simulation <- function(formula=NULL,
                                   baseline.event.rate,
                                   cost.screening=NULL,
                                   cost.keeping=NULL,
                                   estimated.auc,
                                   #roc.type=c("symmetric", "low.tpr.earlier", "high.tpr.earlier"),
                                   roc.type=NULL,
                                   simulation.sample.size=5e+5,
                                   alternative=c("one.sided", "two.sided"),
                                   power=0.9,
                                   alpha=0.025,
                                   reduction.under.treatment=0.3,
                                   n.biomarker.quantiles=20) {

    #############
    ## Check arguments ##
    #############
    if (is.null(baseline.event.rate)) {
        stop("since the data argument was not specified, baseline.event.rate must be specified")
    }
    stopifnot(is.numeric(baseline.event.rate))
    stopifnot(baseline.event.rate > 0 & baseline.event.rate < 1)
    # now we're allowing multiple AUCs to be specified, so we need to do checks element-wise
    if (!(estimated.auc[1] >= 0.5 & estimated.auc[1] <= 1)) {
        stop("Marker 1 must be specified")
    }
    updated.auc <- estimated.auc[!is.na(estimated.auc)]
    n.auc <- length(updated.auc)
    for (i in 1:n.auc) {
        stopifnot(updated.auc[i] >= 0.5 & updated.auc[i] <= 1)
    }
    #roc.type <- match.arg(roc.type)
    updated.roc.type <- roc.type[!is.na(estimated.auc)]
    alternative <- match.arg(alternative)
    if (!all(power > 0 & power < 1)) {
        stop("power should be between 0 and 1")
    }
    if (!all(alpha > 0 & alpha < 1)) {
        stop("alpha should be between 0 and 1")
    }
    if (!all(reduction.under.treatment > 0 & reduction.under.treatment < 1)) {
        stop("reduction.under.treatment should be between 0 and 1")
    }   
    N <- simulation.sample.size
    table.data.list <- vector("list", 3)
    plot.data.check <- NULL
    roc.data.check <- NULL
    for (i in 1:n.auc) {
        simulation.data <- user_auc_and_roc_type_to_data(N=N, baseline.event.rate=baseline.event.rate,
                                                                                auc=updated.auc[i], roc.type=updated.roc.type[i],
                                                                                n.biomarker.quantiles=n.biomarker.quantiles)
        roc.data <- user_auc_to_plots(auc=updated.auc[i], baseline.event.rate=baseline.event.rate, roc.type=updated.roc.type[i], prototypical=FALSE)
        biomarker <- simulation.data$biomarker
        response <- simulation.data$response
        biomarker.quantiles.all <- simulation.data$biomarker.quantiles.all
        selected.biomarker.quantiles <- simulation.data$selected.biomarker.quantiles * 100
###########################################
 ## Calculate the summaries of the vector of biomarker quantiles we want to display ##
###########################################
        if (updated.roc.type[i] %in% c("symmetric", "high.tpr.earlier")) {
            # NNS <- sapply(biomarker.quantiles.all, function(x) N / sum(biomarker > x)) ## 06/21: commenting out since we are no longer using percentiles in controls
            NNS <- 1 / (1 - simulation.data$selected.biomarker.quantiles)
            event.rate <- sapply(biomarker.quantiles.all, function(x) sum(response[biomarker >= x]) / sum(biomarker >= x))
            event.rate[1] <- baseline.event.rate
        } else {
            ## negative biomarker value to keep direction of unequality (since a test-positive is now a lower value of the biomarker)
            # NNS <- sapply(biomarker.quantiles.all, function(x) N / sum(-biomarker > x)) ## 06/21: commenting out since we are no longer using percentiles in controls
            NNS <- 1 / (1 - simulation.data$selected.biomarker.quantiles)
            event.rate <- sapply(biomarker.quantiles.all, function(x) sum(response[-biomarker >= x]) / sum(-biomarker >= x))
            event.rate[1] <- baseline.event.rate
        }
        SS <- sample_size(event.rate=event.rate, reduction.under.treatment=reduction.under.treatment, alpha=alpha, power=power, alternative=alternative)
        N.screen <- SS * NNS # total number of patients needed to be screened
        cost.missing <- is.na(cost.screening) | is.null(cost.screening) | is.na(cost.keeping) | is.null(cost.keeping)
        if (cost.missing == TRUE) {
            plot.data <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.quantiles.all, event.rate, SS, N.screen))
            table.data <- as.data.frame(cbind(paste(selected.biomarker.quantiles, "%", sep=""),
                                              round(event.rate, 2),
                                              round(SS, 0),
                                              round(NNS, 1),
                                              round(N.screen, 0)))
            table.data <- as.data.frame(cbind(apply(table.data[, 1:5], 2, function(x) as.character(x))))
            names(table.data) <- c("Percent of Patients Screened from Trial", "Event Rate Among Biomarker-Positive Patients", "Sample Size", "NNS", "Total Screened")
            print(names(table.data))
            rownames(table.data) <- NULL
        } else {
            print(NNS)
            total.cost <- SS * (cost.keeping + cost.screening * NNS)
            total.cost[1] <- SS[1] * cost.keeping
            cost.reduction.percentage <- ((total.cost[1] - total.cost) / total.cost[1]) * 100
            plot.data <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.quantiles.all, event.rate, SS, N.screen, total.cost, cost.reduction.percentage))
            table.data <- as.data.frame(cbind(paste(selected.biomarker.quantiles , "%", sep=""),
                                              round(event.rate, 2),
                                              round(SS, 0),
                                              round(NNS, 1),
                                              round(N.screen, 0),
                                              round(total.cost, 0),
                                              round(cost.reduction.percentage, 1)))
            table.data <- cbind(apply(table.data[, 1:6], 2, function(x) as.character(x)),
                                paste(as.character(table.data[, 7]), "%", sep=""))
            table.data <- as.data.frame(table.data)
            names(table.data) <- c("Percent of Patients Screened from Trial", "Event Rate Among Biomarker-Positive Patients", "Sample Size", "NNS", "Total Screened", "Total Costs for Screening and Patients in Trial", "Percent Reduction in Total Cost")
            rownames(table.data) <- NULL
        }
        table.data.list[[i]] <- table.data
        #plot.data$AUC <- as.character(updated.auc[i])
        plot.data$Biomarker <- paste("Biomarker", as.character(i), sep=" ")
        plot.data.check <- rbind(plot.data.check, plot.data)
        roc.df <- as.data.frame(cbind(roc.data$fpr.vec, roc.data$tpr.vec))
        roc.df$Biomarker <- paste("Biomarker", as.character(i), sep=" ")
        roc.data.check <- rbind(roc.data.check, roc.df)
    }
    names(roc.data.check) <- c("FPR", "TPR", "Biomarker")
    return(list("plot.data"=plot.data.check,"table.data.list"=table.data.list, "roc.data"=roc.data.check))
}


user_auc_and_roc_type_to_data <- function(N, baseline.event.rate, auc, roc.type, n.biomarker.quantiles) {
    if (roc.type == "low.tpr.earlier") {
        # need to flip case/control labels and eventually call a test "positive" if it is below the threshold, rather than above
        selected.biomarker.quantiles <- seq(from=0, to=0.95, length.out=n.biomarker.quantiles)      
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        biomarker <- numeric(N)
        biomarker[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=1)
        biomarker[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=(1-auc)/auc)
        #biomarker.quantiles.controls <- quantile(-biomarker[response == 0], prob=selected.biomarker.quantiles)  ## commenting out on 06/07/16 to switch to percentiles in all patients
        biomarker.quantiles.all <- quantile(-biomarker, prob=selected.biomarker.quantiles)
    }
    if (roc.type == "symmetric") {
        selected.biomarker.quantiles <- seq(from=0, to=0.95, length.out=n.biomarker.quantiles)
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        sd.cases <- 1
        # solve for (a, b) using Katie's slides 7 and 9, 
        b <- 1 / sd.cases
        a <- sqrt(1 + b^2) * qnorm(auc)
        mean.cases <- a * sd.cases 
        biomarker <- numeric(N)
        biomarker[response == 1] <- rnorm(n=sum(response==1), mean=mean.cases, sd=sd.cases)
        biomarker[response == 0] <- rnorm(n=sum(response==0), mean=0, sd=1)
        #biomarker.quantiles.controls <- quantile(biomarker[response == 0], prob=selected.biomarker.quantiles)
        biomarker.quantiles.all <- quantile(biomarker, prob=selected.biomarker.quantiles)
    } else if (roc.type == "high.tpr.earlier") {
        selected.biomarker.quantiles <- seq(from=0, to=0.95, length.out=n.biomarker.quantiles)
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        biomarker <- numeric(N)
        biomarker[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=1)
        biomarker[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=(1-auc)/auc)
        #biomarker.quantiles.controls <- quantile(biomarker[response == 0], prob=selected.biomarker.quantiles)
        biomarker.quantiles.all <- quantile(biomarker, prob=selected.biomarker.quantiles)
    }
    return(list("selected.biomarker.quantiles"=selected.biomarker.quantiles,
                   "biomarker.quantiles.all"=biomarker.quantiles.all,
                   "response"=response, "biomarker"=biomarker))
}


user.auc.and.roc.type.to.parameters <- function(auc, roc.type) {
    if (roc.type == "low.tpr.earlier") {
        sd.cases <- 2/3
    }
    if (roc.type == "symmetric") {
        sd.cases <- 1
    } else if (roc.type == "high.tpr.earlier") {
        sd.cases <- 3/2
    }
    # solve for (a, b) using Katie's slides 7 and 9, 
    b <- 1 / sd.cases
    a <- sqrt(1 + b^2) * qnorm(auc)
    mean.cases <- a * sd.cases 
    return(list("mean.cases"=mean.cases, "sd.cases"=sd.cases))
}


plot_enrichment_summaries <- function(x,
                                           text.size.x.axis=13,
                                           text.size.y.axis=13,
                                           text.size.plot.title=13,
                                           text.size.axis.ticks=14,
                                           n.col=NULL) {
    plot.data <- x$plot.data
    print(names(plot.data))
    roc.data <- x$roc.data

    cost.indicator <- "total.cost" %in% names(plot.data) & "cost.reduction.percentage" %in% names(plot.data)
    ## ROC curves
    plot.ROC <- ggplot(roc.data, aes(FPR, y=TPR, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) + geom_line(size=0.9) + geom_abline(intercept = 0, slope = 1, linetype="dashed", colour="gray") + coord_fixed() +
        labs(x="FPR", y="TPR") +
        ggtitle("ROC Curve for Specified Biomarkers")  
    plot.ROC <- plot.ROC + expand_limits(y=0) +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title))  + scale_x_continuous(expand = c(0, 0)) +
        theme(legend.text=element_text(size=14), legend.title=element_blank()) +
        theme(legend.key.size = unit(0.80, "cm")) + 
        theme(legend.position=c(0.7, 0.20))
    ## Biomarker percentile vs. sample size
    plot.sample.size <- ggplot(plot.data, aes(selected.biomarker.quantiles, y=SS, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) + geom_line(size=0.9) + geom_point(size=2.5) +
        labs(x="", y="Sample Size") +
        ggtitle("Clinical Trial Total Sample Size") 
    plot.sample.size <- plot.sample.size + expand_limits(y=0) +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title))  + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
    
    ## Biomarker percentile vs. event rate after screening
    plot.event.rate <- ggplot(plot.data, aes(selected.biomarker.quantiles, y=event.rate, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) + geom_line(size=0.9) + geom_point(size=2.5) +
            labs(x="", y="Event Rate") +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title))+
         expand_limits(y=0) + ggtitle("Event Rate Among \n Biomarker-Positive Patients") + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')    
    ## Biomarker percentile vs. total # needing to be screened
    plot.N.screen <- ggplot(plot.data, aes(selected.biomarker.quantiles, y=N.screen, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) 
        if (cost.indicator == TRUE) {
            plot.N.screen <-  plot.N.screen + labs(x="", y="Total Screened")
        } else {
            plot.N.screen <- plot.N.screen + labs(x="Percent of Patients Screened from Trial", y="Total Screened")
        }
    plot.N.screen <- plot.N.screen +
        geom_line(size=0.9) + geom_point(size=2.5) +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title)) +
        ggtitle("Total Number of Patients \n Screened to Enroll Trial") +  scale_x_continuous(expand = c(0, 0)) +
        theme(legend.text=element_text(size=14), legend.title=element_blank()) +
        theme(legend.key.size = unit(1.3, "cm"))
    ## try to plot legend separately (since it's common to all plots)
    one.legend <- g_legend(plot.N.screen)
    plot.N.screen <- plot.N.screen + theme(legend.position = 'none')
    ## Biomarker percentile vs. total cost
    ### Determine whether the total cost information is in the data frame (i.e. whether the user specified cost.screening and cost.retention earlier)
    if (cost.indicator == TRUE) {
        ind.total.cost <- 1:nrow(plot.data)
        plot.total.cost <- ggplot(plot.data[ind.total.cost, ], aes(selected.biomarker.quantiles, y=total.cost, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) + geom_line(size=0.9) + geom_point(size=2.5) + geom_hline(yintercept=0, linetype=2) +
            labs(x="Percent of Patients Screened from Trial", y="Total Cost") +
            ggtitle("Total Costs for Screening \n and Patients in Trial")
        plot.total.cost <- plot.total.cost +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title)) +
            scale_y_continuous(limits = c(0, max(plot.data[, "total.cost"]))) +  scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
        ## Biomarker percentile vs. percentage reduction in total cost (relative to no screening scenario)
        plot.cost.reduction.percentage <- ggplot(plot.data[ind.total.cost, ], aes(selected.biomarker.quantiles, y=cost.reduction.percentage, group=Biomarker, shape=Biomarker, col=Biomarker, linetype=Biomarker, fill=Biomarker)) + geom_line(size=0.9) + geom_point(size=2.5) + geom_hline(yintercept=0, linetype=2) +
        labs(x="Percent of Patients Screened from Trial", y="% Reduction in Total Cost") +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title)) +  scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(limits = c(min(plot.data[, "cost.reduction.percentage"]), max(plot.data[, "cost.reduction.percentage"]))) + ggtitle("Percent Reduction in Total Cost") + theme(legend.position = 'none')
        #args.to.plot <- list("one.legend"=one.legend, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size, "plot.N.screen"=plot.N.screen, "plot.total.cost"=plot.total.cost, "plot.cost.reduction.percentage"=plot.cost.reduction.percentage)
        args.to.plot <- list("plot.ROC"=plot.ROC, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size, "plot.N.screen"=plot.N.screen, "plot.total.cost"=plot.total.cost, "plot.cost.reduction.percentage"=plot.cost.reduction.percentage)
    } else {
        #args.to.plot <- list("one.legend"=one.legend, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size, "plot.N.screen"=plot.N.screen)
        args.to.plot <- list("plot.ROC"=plot.ROC, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size, "plot.N.screen"=plot.N.screen)
    }
    ## Show the plots the user wants to see in a grid (with n.col==2 by default)
    n.plots <- length(names(args.to.plot))
    if (is.null(n.col)) {
        do.call(grid.arrange, c(args.to.plot, list(ncol=floor(n.plots/2))))
    } else {
        do.call(grid.arrange, c(args.to.plot, list(ncol=n.col)))
    }
}

expit <- function(x) {
    result <- exp(x) / (1 + exp(x))
    return(result)
}

user_auc_to_plots <- function(auc, baseline.event.rate, roc.type=NULL, prototypical=TRUE, n=5e+4, n.thresholds=1000, verbose=TRUE) {
    response <- rbinom(n, size=1, prob=baseline.event.rate)
    ## high TPR earlier (orange)
    x.high.tpr.earlier <- rep(NA, n)
    x.high.tpr.earlier[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=1)
    x.high.tpr.earlier[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=(1-auc)/auc)
    result.high.tpr.earlier <- get_roc(x=x.high.tpr.earlier, response=response, test.positive="higher", n.thresholds=n.thresholds, verbose=verbose)
    ## high TPR earlier (blue) 
    x.low.tpr.earlier <- rep(NA, n)
    x.low.tpr.earlier[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=1)
    x.low.tpr.earlier[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=(1-auc)/auc)
    result.low.tpr.earlier <- get_roc(x=x.low.tpr.earlier, response=response, test.positive="lower", n.thresholds=n.thresholds, verbose=verbose)
    ## symmetric (black)
    x.symmetric <- rep(NA, n)
    x.symmetric[response == 0] <- rnorm(n=sum(response==0), mean=0, sd=1)
    x.symmetric[response == 1] <- rnorm(n=sum(response==1), mean=sqrt(2) * qnorm(auc), sd=1)
    result.symmetric <- get_roc(x=x.symmetric, response=response, test.positive="higher", n.thresholds=n.thresholds, verbose=verbose)
    # print estimated AUCs (just for debuggin)
    if (verbose == TRUE) {
        if (prototypical == TRUE) {
            print("this is just the AUC for prototypical ROC curves")
        }
        print(result.high.tpr.earlier$auc.estimate)
        print(result.low.tpr.earlier$auc.estimate)
        print(result.symmetric$auc.estimate)
    }
    # make plots
    if (prototypical == TRUE) {
        par(pty="s")
        plot(result.high.tpr.earlier$fpr.vec, result.high.tpr.earlier$tpr.vec, type="l", lty=2, lwd=3, col="orange", ylim=c(0, 1), xlim=c(0, 1),
             xlab="FPR", ylab="TPR",
             main=paste("Prototypical ROC Curves \n (Example for AUC=", auc, ")", sep=""))
        lines(result.low.tpr.earlier$fpr.vec, result.low.tpr.earlier$tpr.vec, type="l", lty=3, lwd=3, col="cyan")
        lines(result.symmetric$fpr.vec, result.symmetric$tpr.vec, type="l", lty=1, lwd=3, col="black")
        abline(a=0, b=1, lwd=2, lty="dashed", col="gray")
    } else {
        if (roc.type == "symmetric") {
            fpr.vec <- result.symmetric$fpr.vec
            tpr.vec <- result.symmetric$tpr.vec
        } else if (roc.type == "low.tpr.earlier") {
            fpr.vec <- result.low.tpr.earlier$fpr.vec
            tpr.vec <- result.low.tpr.earlier$tpr.vec
        } else if (roc.type == "high.tpr.earlier") {
            fpr.vec <- result.high.tpr.earlier$fpr.vec
            tpr.vec <- result.high.tpr.earlier$tpr.vec
        }
        return(list("fpr.vec"=fpr.vec, "tpr.vec"=tpr.vec))
    }
}

sample_size <- function(event.rate, reduction.under.treatment,
                                 alpha=0.025, power=0.9,
                                 alternative=c("one.sided", "two.sided")) {
    alternative <- match.arg(alternative)
    p1 <- event.rate
    p2 <- event.rate * (1 - reduction.under.treatment)
    z1.one.sided <- qnorm(1 - alpha)
    z1.two.sided <- qnorm(1 - alpha/2)
    z2 <- qnorm(power)
    if (alternative == "one.sided") { 
        SS <- 2 * ( (z1.one.sided * sqrt((p1 + p2) * (1 - (p1 + p2)/2) ) +
                        z2 * sqrt(p1 * (1 - p1) + p2 * (1 - p2)) )^2 / (p1 - p2)^2)
    } else if (alternative == "two.sided") {
        SS <- 2 * ( (z1.two.sided * sqrt((p1 + p2) * (1 - (p1 + p2)/2) ) +
                        z2 * sqrt(p1 * (1 - p1) + p2 * (1 - p2)) )^2 / (p1 - p2)^2)
    } 
    return(SS)
}


# function to compute (FPR, TPR) for a grid of thresholds and estimate AUC
get_roc <- function(x, response, test.positive=c("higher","lower"), n.thresholds=1000, verbose=TRUE) {
    test.positive <- match.arg(test.positive)
    n.healthy <- sum(response == 0)
    n.disease <- sum(response == 1)
    if (test.positive == "lower") {
        x <- -x
    }
    thresholds <- as.numeric(sort(quantile(x=x[response==0], prob= seq(from=0, to=0.99, length.out=n.thresholds), decreasing=FALSE)))
    tpr <- rep(NA, n.thresholds)
    fpr <- rep(NA, n.thresholds)
    x.cases <- x[response==1]
    x.controls <- x[response==0]
    for (i in 1:n.thresholds) {
        tpr[i] <- sum(x.cases > thresholds[i]) / n.disease
        fpr[i] <- sum(x.controls > thresholds[i]) / n.healthy
    }
    if (verbose == TRUE) {
        auc.estimate <- mean(sample(x.cases, size=5e+5, replace=TRUE) > sample(x.controls, size=5e+5, replace=TRUE))
        return(list("auc.estimate"=auc.estimate, "fpr.vec"=fpr, "tpr.vec"=tpr))
    } else {
        return(list("fpr.vec"=fpr, "tpr.vec"=tpr))
    }
}


## trying to plot legend separately
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


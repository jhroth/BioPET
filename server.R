# Jeremy Roth
library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(VGAM)
source("functions.R")



server <- function(input, output) {
    output$roc.choices <- renderPlot({
        user_auc_to_plots(auc=0.75, baseline.event.rate=input$outcome_prevalence)
    })
    summaries <- reactive({
        print(c(as.numeric(input$auc.marker1), as.numeric(input$auc.marker2), as.numeric(input$auc.marker3)))
        summaries <- enrichment_simulation(cost.keeping=input$cost_keeping,
                                           cost.screening=input$cost_screening,
                                           baseline.event.rate=input$outcome_prevalence,
                                           estimated.auc=c(as.numeric(input$auc.marker1), as.numeric(input$auc.marker2), as.numeric(input$auc.marker3)),
                                           roc.type=c(input$roc.type.marker1, input$roc.type.marker2, input$roc.type.marker3),
                                           reduction.under.treatment=input$rate_reduction / 100,
                                           alternative=as.character(input$alternative),
                                           power=as.numeric(input$power),
                                           alpha=as.numeric(input$alpha))
    })
    grid <- function(n.col=2) {
        plot_enrichment_summaries(x=summaries(), n.col=n.col)
    }
    output$summary.table.marker1 <- renderTable({
        print(head(summaries()$table.data.list[[1]]))
        summaries()$table.data.list[[1]]
    })
    output$downloadTableMarker1<- downloadHandler(
        filename = function() {"enrichment_summary_table_biomarker1.csv"},
        content = function(file) {
            my.table <- summaries()$table.data.list[[1]]
            n <- nrow(my.table)
            p <- ncol(my.table)
            print(dim(my.table))
            if (p == 5) {
                my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                     "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve")
                my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                               as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker1), as.character(input$roc.type.marker1))
            } else if (p == 7) {
                my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                     "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve",
                                     "Screening cost", "Retention cost")
                my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                               as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker1), as.character(input$roc.type.marker1),
                               as.character(input$cost_screening), as.character(input$cost_keeping))
            } else{
                stop("invalid inputs specified")
            }
            n.blanks <- n - length(my.inputs.names)
            # rbind(my.inputs, rep("", n.blanks)))
            my.output <- as.data.frame(cbind(my.table, rep("", n), rep("", n), c(my.inputs.names, rep("", n.blanks)), c(my.inputs, rep("", n.blanks))))
            print(names(my.table))
            names(my.output) <- c(names(my.table), "", "", "Input Name", "Input Value")
            write.csv(my.output, file)
        }
     )
    output$summary.table.marker2 <- renderTable({
        my.table.marker2 <- summaries()$table.data.list[[2]]
        my.table.marker2
    })    
    output$downloadTableMarker2<- downloadHandler(
        filename = function() {"enrichment_summary_table_biomarker2.csv"},
        content = function(file) {
            my.table <- summaries()$table.data.list[[2]]
            if (is.null(my.table) == TRUE) {
                write.csv("Biomarker 2 was not specified", file)
            } else {
                n <- nrow(my.table)
                p <- ncol(my.table)
                print(dim(my.table))
                if (p == 6) {
                    my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                         "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve")
                    my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                                   as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker2), as.character(input$roc.type.marker2))
                } else if (p == 8) {
                    my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                         "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve",
                                         "Screening cost", "Retention cost")
                    my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                                   as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker2), as.character(input$roc.type.marker2),
                                   as.character(input$cost_screening), as.character(input$cost_keeping))
                } else{
                    stop("invalid inputs specified")
                }
                n.blanks <- n - length(my.inputs.names)
                                        # rbind(my.inputs, rep("", n.blanks)))
                my.output <- as.data.frame(cbind(my.table, rep("", n), rep("", n), c(my.inputs.names, rep("", n.blanks)), c(my.inputs, rep("", n.blanks))))
                print(names(my.table))
                names(my.output) <- c(names(my.table), "", "", "Input Name", "Input Value")
                write.csv(my.output, file)
            }
        }
    )
    
    output$summary.table.marker3 <- renderTable({
        my.table.marker3 <- summaries()$table.data.list[[3]]
        my.table.marker3
    })    
    output$downloadTableMarker3<- downloadHandler(
        filename = function() {"enrichment_summary_table_biomarker3.csv"},
        content = function(file) {
            my.table <- summaries()$table.data.list[[3]]
            if (is.null(my.table) == TRUE) {
                write.csv("Biomarker 3 was not specified", file)
            } else {
                n <- nrow(my.table)
                p <- ncol(my.table)
                print(dim(my.table))
                if (p == 6) {
                    my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                         "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve")
                    my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                                   as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker3), as.character(input$roc.type.marker3))
                } else if (p == 8) {
                    my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                         "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve",
                                         "Screening cost", "Retention cost")
                    my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                                   as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker3), as.character(input$roc.type.marker3),
                                   as.character(input$cost_screening), as.character(input$cost_keeping))
                } else{
                    stop("invalid inputs specified")
                }
                n.blanks <- n - length(my.inputs.names)
                                        # rbind(my.inputs, rep("", n.blanks)))
                my.output <- as.data.frame(cbind(my.table, rep("", n), rep("", n), c(my.inputs.names, rep("", n.blanks)), c(my.inputs, rep("", n.blanks))))
                print(names(my.table))
                names(my.output) <- c(names(my.table), "", "", "Input Name", "Input Value")
                write.csv(my.output, file)
            }
        }
    )
    
    output$grid <- renderPlot({
        grid()
    })
    output$downloadPlots <- downloadHandler(
        filename = function() {"enrichment_summary_plots.png"},
        content = function(file) {
            png(file, width=720, height=720)
            grid(n.col=2)
            dev.off()
        }
     )
}

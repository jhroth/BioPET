# Jeremy Roth
#rm(list=ls())
library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(VGAM)
source("functions.R")

ui <- fluidPage(
    h1("BioPET:  Biomarker Prognostic Enrichment Tool"),
    tabsetPanel(
        tabPanel("About BioPET",
             fluidRow(
                 column(width=12, offset=0,
                        HTML(" <br> 
<p> <h4> <strong>BioPET </strong> is a tool for helping investigators evaluate  whether a biomarker or risk model is useful for prognostic enrichment of a clinical trial.  The clinical trial will study the efficacy of an intervention for reducing the rate of an unwanted event.  As described by Temple (2010), a biomarker or risk model can allow the intervention to be evaluated on a study population with higher rates of the event than the general population.  This decreases the necessary sample size for the clinical trial, which can have both practical and ethical advantages. </h4> </p>
<p> <h4> <strong> BioPET </strong> is also available as a package for the R Statistical Computing Platform.  The R package offers extended functionality compared to the webtool.  In particular, the R package allows investigators to analyze their own biomarker data rather than relying on prototypical ROC curves. </h4> </p> <br>
<p> <h5> <strong> Reference </strong> </h5>
<p> Temple R. Enrichment of clinical study populations. Clin Pharmacol Ther 2010; 88: 774-778. </h4> </p>
<br> <hr>")
                 )
            ),
            fluidRow(
                column(width=9, offset=0,
                            h2(strong("Input Information"), align="center")
                          ),
                column(width=3, offset=0,
                            h2(strong("Results"), align="center")
                          )),
             fluidRow(
                 column(width=3, offset=0,
                        h3(strong("Clinical Trial Information"), align="center"),
                        br(),
                        HTML("<p> <h5> <strong> Event rate in the non-intervention group </strong> </p> <p> The prevalence of the adverse event in the non-intervention group. For example, 20% of patients in the non-intervention group may be expected to experience the unwanted event. </h5> </p> <br>
                                  <p> <h5> <strong> Percent reduction in event rate under treatment </strong> </p> <p> The effect size of the intervention as represented by the percent reduction in the event rate for patients using the therapy that the trial should be powered to detect. For example, we may want to design a trial powered to detect a 30% reduction in the event rate. </h5> </p> <br>
                                  <p> <h5> <strong> Form of alternative hypothesis </strong> </p> <p> Indicates whether the study will use one- or two-sided hypothesis testing. </h5> </p> <br>
                                  <p> <h5> <strong> Type I error rate </strong> </p> <p> The probability of rejecting the null hypothesis, given that the null hypothesis is actually true. Common settings are 0.025 and 0.05.</h5> </p> <br>
                                  <p> <h5> <strong> Power </strong> </p> <p> The probability of rejecting the null hypothesis, given that the null hypothesis is actually false. For example, we might design our clinical trial to have 90% power to detect the treatment effect. </h5> </p>")
                        ),
                 column(width=3, offset=0,
                        h3(strong("Biomarker Information"), align="center"),
                        br(),
                        HTML("<p> <h5> <strong> AUC of biomarker </strong> </p> <p> The area under the ROC curve for the biomarker, summarizing the biomarker's ability to distinguish between cases and controls. </h5> </p> <br>
                             <p> <h5> <strong> ROC curve for biomarker </strong> </p> <p> Asks for additional information about the shape of the ROC curve that yields the inputted AUC. The three displayed ROC curves have the same AUC, but have different shapes (symmetric, left-shifted, or right-shifted) which provide information about the predictive capacility of the biomarker. </h5> </p>")
                 ),
                 column(width=3, offset=0,
                        h3(strong("Cost Information"), align="center"),
                        h4("(optional)", align="center"),
                        HTML("<p> <h5> <strong> Cost of screening a patient to determine trial eligibility </strong> </p> <p> For example the cost of measuring the biomarker to determine patient eligibility for the trial may be $100. </h5> </p> <br>
                                    <p> <h5> <strong> “Cost of running a patient through the trial”  </strong> </p> <p> For example, the cost of enrolling and retaining a patient in a trial may be $1000.</h5> </p>")
                 ),
                 column(width=3, offset=0,
                        br(),
                        br(),
                        br(),
                        br(),
                        HTML("<p> <h5> <strong> Sample size </strong> </p> <p> The sample size required for a clinical trial enrolling only patients who are biomarker-positive. </h5> </p> <br>
                                  <p> <h5> <strong> Number needed to screen (NNS) </strong> </p> <p> The estimated number of patients who need to be screened to identify one patient eligible for the trial. </h5> </p> <br>
                                   <p> <h5> <strong> “Event rate among biomarker-positive patients </strong> </p> <p> The estimated event rate among the trial participants if the biomarker were used for prognostic enrichment. </h5> </p> <br>
                                   <p> <h5> <strong> Total screened </strong> </p> <p> The estimated total number of individuals who must be screened to enroll the prognostically enriched trial. </h5> </p> <br>
                                   <p> <h5> <strong> Total cost </strong> </p> <p>  The estimated total cost of running the trial if the biomarker were used for prognostic enrichment. </h5> </p>")
                        )
             )
        ),
    tabPanel("Inputs", 
    fluidRow(
         column(width=3,
               h3(strong("Clinical Trial Information"), align="center"),
               br(),
               numericInput(inputId="outcome_prevalence",
                            label="Event rate in the non-intervention group",
                            value=0.20, min=0.01, max=0.99, step=0.01),
               numericInput(inputId="rate_reduction",
                            label="Percent reduction in event rate under treatment",
                            value=30, min=1, max =99, step=1),
               radioButtons(inputId="alternative",
                            label=HTML("Form of alternative hypothesis"),
                            choices=c("one-sided"="one.sided", "two-sided"="two.sided")),
               radioButtons(inputId="alpha",
                            label=HTML("Type I error rate (&alpha;)"),
                            choices=c(0.01, 0.025, 0.05),
                            selected=0.025),
               radioButtons(inputId="power",
                            label=HTML("Power (1-&beta;)"),
                            choices=c("90%"=0.90, "80%"=0.80, "70%"=0.70))                       
        ),
        column(width=6,
               h3(strong("Biomarker Information"), align="center"),
               column(width=4, 
                      numericInput(inputId="auc.marker1",
                                   label=h5(strong("AUC of Biomarker 1 \n (required)"), align="center"),
                                   value=0.7, min=0.51, max=0.99, step=0.01),
                      radioButtons(inputId="roc.type.marker1",
                                   label="Which of the curves below looks most like Biomarker 1's ROC curve?",
                                   choices=c("left-shifted (orange, dashed curve)"="high.tpr.earlier",
                                             "symmetric (black, solid curve)"="symmetric",
                                             "right-shifted (cyan, dotted curve)"="low.tpr.earlier"),
                                   selected="symmetric")
                      ),
               column(width=4,
                      numericInput(inputId="auc.marker2",
                                   label=h5(strong("AUC of Biomarker 2 \n (optional)"), align="center"),
                                   value=NA, min=0.51, max=0.99, step=0.01),
                      radioButtons(inputId="roc.type.marker2",
                                   label="Which of the curves below looks most like Biomarker 2's ROC curve?",
                                   choices=c("left-shifted (orange, dashed curve)"="high.tpr.earlier", "symmetric (black, solid curve)"="symmetric", "right-shifted (cyan, dotted curve)"="low.tpr.earlier"),
                                      selected="symmetric")
#                                   selected=character(0))
                      ),
               column(width=4, 
                      numericInput(inputId="auc.marker3",
                                   label=h5(strong("AUC of Biomarker 3 \n (optional)"), align="center"),
                                   value=NA, min=0.51, max=0.99, step=0.01),
                      radioButtons(inputId="roc.type.marker3",
                                   label="Which of the curves below looks most like Biomarker 3's ROC curve?",
                                   choices=c("left-shifted (orange, dashed curve)"="high.tpr.earlier", "symmetric (black, solid curve)"="symmetric", "right-shifted (cyan, dotted curve)"="low.tpr.earlier"),
                                   selected="symmetric")
#                                   selected=character(0))
                      ),
               plotOutput("roc.choices")
        ),
        column(width=3,
               h3(strong("Cost Information"), align="center"),
               h4("(optional)", align="center"),
               numericInput(inputId="cost_screening",
                            label="Cost of screening a patient to determine trial eligiblity",
                            value=NULL, min=1, max=1e+05, step=0.5),
               numericInput(inputId="cost_keeping",
                            label="Cost of running a patient through the trial",
                            value=NULL, min=1, max=1e+05, step=0.5),
               br(),
               br(),
               submitButton("Submit", icon("refresh"), width='100%')
               ))),
            tabPanel("Results Summary",
                      tabsetPanel(
                          tabPanel("Biomarker 1",
                                  fluidRow(
                                      column(width=8,
                                             h4(strong("Table 1: Summary Measures as a Function of Biomarker 1 Value Used for Screening")),
                                             div(tableOutput("summary.table.marker1"), style="font-size:100%")
                                             ),
                                      column(width=4,
                                             br(),
                                             downloadButton(outputId="downloadTableMarker1", label="Download results in a .csv file")
                                             )
                                  )
                              ),
                          tabPanel("Biomarker 2",
                                  fluidRow(
                                      column(width=8,
                                             h4(strong("Table 2: Summary Measures as a Function of Biomarker 2 Value Used for Screening")),
                                             div(tableOutput("summary.table.marker2"), style="font-size:100%")
                                             ),
                                      column(width=4,
                                             br(),
                                             downloadButton(outputId="downloadTableMarker2", label="Download results in a .csv file")
                                             )
                                  )
                              ),
                          tabPanel("Biomarker 3",
                                  fluidRow(
                                      column(width=8,
                                             h4(strong("Table 3: Summary Measures as a Function of Biomarker 3 Value Used for Screening")),
                                             div(tableOutput("summary.table.marker3"), style="font-size:100%")
                                             ),
                                      column(width=4,
                                             br(),
                                             downloadButton(outputId="downloadTableMarker3", label="Download results in a .csv file")
                                             )
                                  )
                              )
                          )
                     ),
            tabPanel("Plots",
                 fluidRow(
                     column(width=8, offset=0,
                            #plotOutput("grid", width = "100%", height="600px"))
                            plotOutput("grid", width="90%", height="800px")
                     ),
                     column(width=4, offset=0,
                            br(),
                            downloadButton(outputId="downloadPlots", label="Download plots in a .png file")
                    )
               )
        )
    )
)


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
            if (p == 6) {
                my.inputs.names <- c("Baseline event rate", "Percent reduction in event rate under treatment",
                                     "Form of alternative hypothesis", "Type I error rate", "Power", "AUC", "Biomarker ROC curve")
                my.inputs <- c(as.character(input$outcome_prevalence), as.character(input$rate_reduction),
                               as.character(input$alternative), as.character(input$alpha), as.character(input$power), as.character(input$auc.marker1), as.character(input$roc.type.marker1))
            } else if (p == 8) {
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

shinyApp(ui = ui, server = server)

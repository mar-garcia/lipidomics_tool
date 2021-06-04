library(shiny)
library(MetaboCoreUtils)
library(OrgMassSpecR)

.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

sn <- data.frame(
  "C" = rep(seq(14, 22), each = 4),
  "db" = rep(seq(0, 3), 9)
)
sn$formula <- paste0("C", sn$C, "H", sn$C*2 - 2*sn$db, "O2")
sn$sn <- paste0(sn$C, ":", sn$db)
sn$mass <- NA
sn.list <- list()
for(i in seq(nrow(sn))){
  sn$mass[i] <- MonoisotopicMass(formula = ListFormula(sn$formula[i]))
  
  sn.list[[i]] <- sn$sn[i]
  names(sn.list)[[i]] <- sn$sn[i]
}

# UI ------------------------------------------------------------------------
ui <- navbarPage(
  "Lipidomics",
  
  ## PA ----
  tabPanel(
    "Phosphatidic acids (PAs)",
    column(4, h3("Formula"),
           fluidRow(
             column(2, numericInput("paC", "C", value = 32)),
             column(2, numericInput("padb", "db", value = 0))
           ),
           column(4, fluidRow(verbatimTextOutput("paformula"))),
           fluidRow(),
           fluidRow(h3("m/z values"), verbatimTextOutput("pamzvals"))
    ),
    column(1), 
    column(6, h3("MS2 (-)"),
           fluidRow(
             column(3, numericInput("paion1", "ion1", value = 0)),
             column(3, numericInput("paion2", "ion2", value = 0))
           ),
           fluidRow(
             column(3, verbatimTextOutput("pasn1")),
             column(3, verbatimTextOutput("pasn2")),
             column(3, verbatimTextOutput("pasum"))
           ),
           hr(),
           fluidRow(
             column(4, selectInput("paion1x", "sn1", choices = sn.list)),
             column(4, selectInput("paion2x", "sn2", choices = sn.list))
           ),
           fluidRow(
             column(4, verbatimTextOutput("pasn1x")),
             column(4, verbatimTextOutput("pasn2x"))
           )
    )
  ), # close tab PAs
  
  
  ## PA methylated ----
  tabPanel(
    "Methylated phosphatidic acids (mPAs)",
    column(4, h3("Formula"),
           fluidRow(
             column(2, numericInput("mpaC", "C", value = 32)),
             column(2, numericInput("mpadb", "db", value = 0))
           ),
           column(4, fluidRow(verbatimTextOutput("mpaformula"))),
           fluidRow(),
           fluidRow(h3("m/z values"), verbatimTextOutput("mpamzvals"))
    ),
    column(1), 
    column(6, h3("MS2 (-)"),
           fluidRow(
             column(3, numericInput("mpaion1", "ion1", value = 0)),
             column(3, numericInput("mpaion2", "ion2", value = 0))
           ),
           fluidRow(
             column(3, verbatimTextOutput("mpasn1")),
             column(3, verbatimTextOutput("mpasn2"))
           )),
    fluidRow(),
    hr(),
    h3("Commonly occuring product ions for mPAs:"),
    column(3,
           strong("Positive [M+NH4]+:"),
           tags$li("[M + H]+"),
           tags$li("[M + H - methyl-phosphate]+")),
    column(3, 
           strong("Negative [M-H]-:"),
           tags$li("[sn1 - H]-"),
           tags$li("[sn2 - H]-"),
           tags$li("[M - H - (sn1-H2O)]-"))
  ), # close mPA
  
  
  ## DAG ----
  tabPanel(
    "Diacylglycerols (DAGs)",
    column(4, h3("Formula"),
           fluidRow(
             column(2, numericInput("dagC", "C", value = 39)),
             column(2, numericInput("dagdb", "db", value = 0))
           ),
           column(4, fluidRow(verbatimTextOutput("dagformula"))),
           fluidRow(),
           fluidRow(h3("m/z values"), verbatimTextOutput("dagmzvalspos")),
           fluidRow("", verbatimTextOutput("dagmzvalsneg"))), 
    column(1), 
    column(6, h3("MS2 (+)"),
           fluidRow(
             column(3, numericInput("dagion1", "ion1", value = 0)),
             column(3, numericInput("dagion2", "ion2", value = 0))
           ),
           fluidRow(
             column(3, verbatimTextOutput("dagsn1")),
             column(3, verbatimTextOutput("dagsn2")),
             column(3, verbatimTextOutput("dagsum"))
           ))), # close tab DAG
  
  ## TAG ----
  tabPanel(
    "Triacylglycerols (TAGs)",
    column(2, h3("Formula"),
           fluidRow(
             column(6, numericInput("tagC", "C", value = 54)),
             column(6, numericInput("tagdb", "db", value = 0))
           ),
           fluidRow(verbatimTextOutput("tagformula")),
           fluidRow(h3("m/z values"), verbatimTextOutput("tagmzvals"))),
    column(1),
    column(6, h3("MS2 (+)"),
           fluidRow(
             column(3, numericInput("tagion1", "ion1", value = 0)),
             column(3, numericInput("tagion2", "ion2", value = 0)),
             column(3, numericInput("tagion3", "ion3", value = 0))
           ),
           fluidRow(
             column(3, verbatimTextOutput("tagsn1")),
             column(3, verbatimTextOutput("tagsn2")),
             column(3, verbatimTextOutput("tagsn3")),
             column(3, verbatimTextOutput("tagsum"))
           ),
           hr(),
           fluidRow(
             column(4, selectInput("tagsn1", "sn1", choices = sn.list)),
             column(4, selectInput("tagsn2", "sn2", choices = sn.list)),
             column(4, selectInput("tagsn3", "sn3", choices = sn.list))
           ),
           fluidRow(fluidRow(verbatimTextOutput("tagms2")))
    )
    
  ) # close tab TAG
)# close ui

# SERVER ---------------------------------------------------------------------
server <- function(input, output) {
  
  ## PA ----
  pafml <- reactive({
    paste0("C", input$paC + 3, "H", 
           input$paC*2 - (2 + 2*input$padb) + 7, "O8P")
  })
  
  pamass <- reactive({
    MonoisotopicMass(formula = ListFormula(pafml()))
  })
  
  output$paformula <- renderPrint({pafml()})
  
  output$pamzvals <- renderPrint({
    tmp <- unlist(mass2mz(
      pamass(), 
      adduct = c("[M+H]+", "[M+NH4]+", "[2M+H]+", "[2M+NH4]+", "[M-H]-")))
    tmp2 <- colnames(tmp)
    tmp <- c(as.numeric(unlist(mass2mz(pamass(), "[M+H]+"))) - 
               MonoisotopicMass(formula = ListFormula("H3PO4")), tmp)
    names(tmp) <- c("[M+H-PA]+", tmp2)
    tmp
  })
  
  output$pasn1 <- renderPrint({
    paste0(
      sn$sn[unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion1, 
        sn$mass, ppm = 10))], " (",
      round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - 
              sn$mass[unlist(matchWithPpm(
                unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion1, 
                sn$mass, ppm = 10))] + 
              MonoisotopicMass(formula = ListFormula("H2O")), 4), ")"
    )
  })
  
  output$pasn2 <- renderPrint({
    paste0(
      sn$sn[unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion2, 
        sn$mass, ppm = 10))], " (",
      round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - 
              sn$mass[unlist(matchWithPpm(
                unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion2, 
                sn$mass, ppm = 10))] + 
              MonoisotopicMass(formula = ListFormula("H2O")), 4), ")"
    )
  })
  
  output$pasum <- renderPrint({
    paste0(
      sn$C[unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion1, 
        sn$mass, ppm = 10))] +
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion2, 
          sn$mass, ppm = 10))],
      ":",
      sn$db[unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion1, 
        sn$mass, ppm = 10))] +
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion2, 
          sn$mass, ppm = 10))]
    )
  })
  
  output$pasn1x <- reactive({
    round(unlist(mass2mz(sn$mass[sn$sn == input$paion1x], "[M-H]-")), 5)
  })
  
  output$pasn2x <- reactive({
    round(unlist(mass2mz(sn$mass[sn$sn == input$paion2x], "[M-H]-")), 5)
  })
  
  
  ## PA methylated -----
  
  mpafml <- reactive({
    paste0("C", input$mpaC + 4, "H", 
           input$mpaC*2 - (2 + 2*input$mpadb) + 9, "O8P")
  })
  
  mpamass <- reactive({
    MonoisotopicMass(formula = ListFormula(mpafml()))
  })
  
  output$mpaformula <- renderPrint({mpafml()})
  
  output$mpamzvals <- renderPrint({
    tmp <- unlist(mass2mz(
      mpamass(), 
      adduct = c("[M+H]+", "[M+NH4]+", "[2M+H]+", "[2M+NH4]+", "[M-H]-", "[2M-H]-")))
    tmp2 <- colnames(tmp)
    tmp <- c(as.numeric(unlist(mass2mz(mpamass(), "[M+H]+"))) - 
               MonoisotopicMass(formula = ListFormula("H3PO4CH2")), tmp)
    names(tmp) <- c("[M+H-mPA]+", tmp2)
    tmp
  })
  
  output$mpasn1 <- renderPrint({
    paste0(
      sn$sn[unlist(matchWithPpm(input$mpaion1 + 1.007276, sn$mass, ppm = 10))], " (",
      round(
        unlist(mass2mz(mpamass(), adduct = c("[M-H]-"))) - 
          MonoisotopicMass(formula = ListFormula(subtractElements(
            sn$formula[unlist(matchWithPpm(input$mpaion1 + 1.007276, 
                                           sn$mass, ppm = 10))], "H2O"))), 
        4), ")"
    )
  })
  
  output$mpasn2 <- renderPrint({
    paste0(
      sn$sn[unlist(matchWithPpm(input$mpaion2 + 1.007276, sn$mass, ppm = 10))], " (",
      round(
        unlist(mass2mz(mpamass(), adduct = c("[M-H]-"))) - 
          MonoisotopicMass(formula = ListFormula(subtractElements(
            sn$formula[unlist(matchWithPpm(input$mpaion2 + 1.007276, 
                                           sn$mass, ppm = 10))], "H2O"))), 
        4), ")"
    )
  })
  
  
  ## DAG ----
  
  dagfml <- reactive({
    paste0("C", input$dagC + 3, "H", 
           input$dagC*2 - (2 + 2*input$dagdb) + 6, "O5")
  })
  
  dagmass <- reactive({
    MonoisotopicMass(formula = ListFormula(dagfml()))
  })
  
  output$dagformula <- renderPrint({dagfml()})
  
  output$dagmzvalspos <- renderPrint({
    tmp <- unlist(mass2mz(
      dagmass(), 
      adduct = c("[M+H-H2O]+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[2M+NH4]+", "[2M+Na]+")))
    tmp2 <- colnames(tmp)
    tmp <- c(tmp, as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
               MonoisotopicMass(formula = ListFormula("C2H7N")))
    names(tmp) <- c(tmp2, "[M+C2H8N]+")
    tmp <- tmp[c(1:4, 7, 5:6)]
    tmp
  })
  
  output$dagmzvalsneg <- renderPrint({
    unlist(mass2mz(
      dagmass(), 
      adduct = c("[M+CHO2]-", "[2M+CHO2]-")))
  })
  
  output$dagms2 <- renderPrint({
    ms2()
  })
  
  output$dagsn1 <- renderPrint({
    sn$sn[unlist(matchWithPpm(
      unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion1, 
      sn$mass, ppm = 10))]
  })
  
  output$dagsn2 <- renderPrint({
    sn$sn[unlist(matchWithPpm(
      unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion2, 
      sn$mass, ppm = 10))]
  })
  
  output$dagsum <- renderPrint({
    paste0(
      sn$C[unlist(matchWithPpm(
        unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion1, 
        sn$mass, ppm = 10))] +
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion2, 
          sn$mass, ppm = 10))],
      ":",
      sn$db[unlist(matchWithPpm(
        unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion1, 
        sn$mass, ppm = 10))] +
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(dagmass(), adduct = c("[M+H]+"))) - input$dagion2, 
          sn$mass, ppm = 10))]
    )
  })
  
  ## TAG ----
  
  tagfml <- reactive({
    paste0("C", input$tagC + 3, "H", 
           input$tagC*2 - (3 + 2*input$tagdb) + 5, "O6")
  })
  
  tagmass <- reactive({
    MonoisotopicMass(formula = ListFormula(tagfml()))
  })
  
  tagms2 <- reactive({
    tmp <- c(
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$sn1],
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$sn2],
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$sn3])
    tmp <- unique(tmp)
    names(tmp) <- unique(paste0("[M+H-C", c(input$sn1, input$sn2, input$sn3), "]+"))
    tmp
  })
  
  output$tagformula <- renderPrint({tagfml()})
  
  output$tagmzvals <- renderPrint({
    unlist(mass2mz(tagmass(), adduct = c("[M+H]+", "[M+NH4]+")))
  })
  
  output$tagms2 <- renderPrint({
    tagms2()
  })
  
  output$tagsn1 <- renderPrint({
    sn$sn[unlist(matchWithPpm(
      unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion1, 
      sn$mass, ppm = 10))]
  })
  
  output$tagsn2 <- renderPrint({
    sn$sn[unlist(matchWithPpm(
      unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion2, 
      sn$mass, ppm = 10))]
  })
  
  output$tagsn3 <- renderPrint({
    sn$sn[unlist(matchWithPpm(
      unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion3, 
      sn$mass, ppm = 10))]
  })
  
  output$tagsum <- renderPrint({
    paste0(
      sn$C[unlist(matchWithPpm(
        unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion1, 
        sn$mass, ppm = 10))] +
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion2, 
          sn$mass, ppm = 10))] +
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion3, 
          sn$mass, ppm = 10))],
      ":",
      sn$db[unlist(matchWithPpm(
        unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion1, 
        sn$mass, ppm = 10))] +
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion2, 
          sn$mass, ppm = 10))] +
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(tagmass(), adduct = c("[M+H]+"))) - input$tagion3, 
          sn$mass, ppm = 10))]
    )
  })
  
} # close server

shinyApp(ui = ui, server = server)
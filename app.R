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

mzdif.pos <- data.frame(rbind(
  c(MonoisotopicMass(formula = ListFormula("NH3")), "loss NH3 -> PA / mPA / PI / DAG / TAG / MGDG / DGDG"),
  c(MonoisotopicMass(formula = ListFormula("H2O")), "loss H2O -> Lyso PC"),
  c(MonoisotopicMass(formula = ListFormula("NH3H2O")), "loss NH3 & H2O -> DAG"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4NH3")), "loss NH3 & phosphate -> PA"),
  c(MonoisotopicMass(formula = ListFormula("NH3H3PO4CH2")), "loss NH3 & mPA -> mPA"),
  c(MonoisotopicMass(formula = ListFormula("C3H9N")), "loss C3H9N -> PC / Carnitine"),
  c(MonoisotopicMass(formula = ListFormula("C5H14NO4P")), "loss C5H14NO4P -> PC"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) - 0.984, 
    "loss phosphoinositol (NH4 ion) -> PI"),
  c(MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 
    "loss phosphoethanolamine -> PE"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> MGDG [M+Na]+")),
  cbind(sn$mass + MonoisotopicMass(formula = ListFormula("NH3")), 
        paste("loss NH3 &", sn$sn, "-> DAG / TAG")),
  
  c(MonoisotopicMass(formula = ListFormula("C6H12O6")) - 0.984, 
    "loss 1 galactose (NH4 ion) -> MGDG"),
  c(MonoisotopicMass(formula = ListFormula("C6H12O6")) + 
      MonoisotopicMass(formula = ListFormula("NH3")), 
    "loss NH4 & 1 galactose -> MGDG"),
  c(MonoisotopicMass(formula = ListFormula("C12H22O11")) - 0.984, 
    "loss 2 galactoses (NH4 ion) -> DGDG"),
  c(MonoisotopicMass(formula = ListFormula("C12H22O11")) + 
      MonoisotopicMass(formula = ListFormula("NH3")), 
    "loss NH4 & 2 galactoses -> DGDG"),
  cbind(MonoisotopicMass(formula = ListFormula("C5H15NO4P")), "[Phosphocholine]+ -> Lyso PC")
))
colnames(mzdif.pos) <- c("dif", "add")
mzdif.pos$dif <- as.numeric(mzdif.pos$dif)

mzdif.neg <- data.frame(rbind(
  c(MonoisotopicMass(formula = ListFormula("HCOOH")), "loss HCOOH -> MGDG / DGDG"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> PA / PG / PI")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), 
        paste0("loss '", sn$sn, "-H2O' -> PA / mPA / PG / PE")),
  cbind((sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))) + 
          MonoisotopicMass(formula = ListFormula("C3H8O3")), 
        paste0("loss '", sn$sn, "-H2O' & glycerol -> PA / PG / PE")),
  cbind(sn$mass + MonoisotopicMass(formula = ListFormula("HCOOH")), 
        paste("loss HCOOH &", sn$sn, "-> DGDG")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))
        + MonoisotopicMass(formula = ListFormula("HCOOH")), 
        paste0("loss HCOOH & '", sn$sn, "-H2O' -> MGDG")),
  cbind(MonoisotopicMass(formula = ListFormula("HCOOHCH2")), 
        "loss HCOOH & CH2 -> Lyso PC / PC"),
  c(MonoisotopicMass(formula = ListFormula("HCOOHC7H13NO2")), "loss HCOOH & C7H13NO2 -> Carnitine"),
  cbind(mass2mz(sn$mass, "[M-H]-"), paste0("[", sn$sn, "-H]- -> PA / mPA / PG / PE / Lyso PC / PI / MGDG"))
))
colnames(mzdif.neg) <- c("dif", "add")
mzdif.neg$dif <- as.numeric(mzdif.neg$dif)

# UI ------------------------------------------------------------------------
ui <- navbarPage(
  "Lipidomics",
  
  tabPanel(
    "MS2",
    column(5, h2("ESI+"),
           numericInput("posprec", "Precursor", value = 0),
           fluidRow(
             column(2, numericInput("posfrag1", "Fragment 1", value = 0)),
             column(8, verbatimTextOutput("posfrag1add"))
           ),
           fluidRow(
             column(2, numericInput("posfrag2", "Fragment 2", value = 0)),
             column(8, verbatimTextOutput("posfrag2add"))
           ),
           fluidRow(
             column(2, numericInput("posfrag3", "Fragment 3", value = 0)),
             column(8, verbatimTextOutput("posfrag3add"))
           ),
           fluidRow(
             column(2, numericInput("posfrag4", "Fragment 4", value = 0)),
             column(8, verbatimTextOutput("posfrag4add"))
           )),
    column(1),
    column(5, h2("ESI-"),
           numericInput("negprec", "Precursor", value = 0),
           fluidRow(
             column(2, numericInput("negfrag1", "Fragment 1", value = 0)),
             column(8, verbatimTextOutput("negfrag1add"))
           ),
           fluidRow(
             column(2, numericInput("negfrag2", "Fragment 2", value = 0)),
             column(8, verbatimTextOutput("negfrag2add"))
           ),
           fluidRow(
             column(2, numericInput("negfrag3", "Fragment 3", value = 0)),
             column(8, verbatimTextOutput("negfrag3add"))
           ),
           fluidRow(
             column(2, numericInput("negfrag4", "Fragment 4", value = 0)),
             column(8, verbatimTextOutput("negfrag4add"))
           ),
           fluidRow(
             column(2, numericInput("negfrag5", "Fragment 5", value = 0)),
             column(8, verbatimTextOutput("negfrag5add"))
           ),
           fluidRow(
             column(2, numericInput("negfrag6", "Fragment 6", value = 0)),
             column(8, verbatimTextOutput("negfrag6add"))
           )
    ),
  ), # close "MS2" tab
  
  ## GP ----
  navbarMenu(
    "Glycerophospholipids (GPs)",
    ### PA ----
    tabPanel(
      "Phosphatidic acids (PAs)",
      h1("Phosphatidic acids (PAs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("paC", "C", value = 32)),
               column(2, numericInput("padb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("paformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("pamzvals1")),
             fluidRow(verbatimTextOutput("pamzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("pafragpos")),
             fluidRow(h4("ESI-"), 
                      fluidRow(
                        column(3, numericInput("paion1", "ion1", value = 0)),
                        column(3, numericInput("paion2", "ion2", value = 0))
                      ),
                      fluidRow(
                        column(3, htmlOutput("pasn1")),
                        column(3, htmlOutput("pasn2")),
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
             )),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for PAs:"),
      column(3,
             strong("Positive [M+NH4]+:"),
             tags$li("[M + H]+"),
             tags$li("[M + H - phosphate]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[M - H - sn1]-"),
             tags$li("[M - H - (sn1-H2O)]-"),
             tags$li("[sn1 - H]-"),
             tags$li("[M - H - sn2]-"),
             tags$li("[M - H - (sn2-H2O)]-"),
             tags$li("[sn2 - H]-")
             
      )
    ), # close tab PAs
    
    
    ### PA methylated ----
    tabPanel(
      "Methylated phosphatidic acids (mPAs)",
      h1("Methylated phosphatidic acids (mPAs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("mpaC", "C", value = 32)),
               column(2, numericInput("mpadb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("mpaformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("mpamzvals1")),
             fluidRow(verbatimTextOutput("mpamzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("mpafragpos")),
             fluidRow(h4("ESI-"), 
                      fluidRow(
                        column(3, numericInput("mpaion1", "ion1", value = 0)),
                        column(3, numericInput("mpaion2", "ion2", value = 0))
                      ),
                      fluidRow(
                        column(3, verbatimTextOutput("mpasn1")),
                        column(3, verbatimTextOutput("mpasn2")),
                        column(3, verbatimTextOutput("mpasum"))
                      ))),
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
      
      ### PG ----
      tabPanel(
        "Phosphatidylglycerols (PGs)",
        h1("Phosphatidylglycerols (PGs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("pgC", "C", value = 32)),
                 column(2, numericInput("pgdb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("pgformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("pgmzvals1")),
               fluidRow(verbatimTextOutput("pgmzvals2"))
        ),
        column(1), 
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("pgfragpos")),
               fluidRow(h4("ESI-"), 
                        fluidRow(
                          column(3, numericInput("pgion1", "ion1", value = 0)),
                          column(3, numericInput("pgion2", "ion2", value = 0))
                        ),
                        fluidRow(
                          column(3, htmlOutput("pgsn1")),
                          column(3, htmlOutput("pgsn2")),
                          column(3, verbatimTextOutput("pgsum"))
                        )
               )
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for PGs:"),
        column(3,
               strong("Positive:")),
        column(3, 
               strong("Negative [M-H]-:"),
               tags$li("[M - H - (sn1-H2O)]-"),
               tags$li("[M - H - sn1]-"),
               tags$li("[M - H - (sn1-H2O) - glycerol]-"),
               tags$li("[sn1 - H]-"),
               br(),
               tags$li("[M - H - (sn2-H2O)]-"),
               tags$li("[M - H - sn2]-"),
               tags$li("[M - H - (sn2-H2O) - glycerol]-"),
               tags$li("[sn2 - H]-")
               
        )
      ), # close tab PGs
      
      ### PE ----
      tabPanel(
        "Phosphatidylethanolamines (PEs)",
        h1("Phosphatidylethanolamines (PEs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("peC", "C", value = 32)),
                 column(2, numericInput("pedb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("peformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("pemzvals1")),
               fluidRow(verbatimTextOutput("pemzvals2"))
        ),
        column(1), 
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("pefragpos")),
               fluidRow(h4("ESI-"), 
                        fluidRow(
                          column(3, numericInput("peion1", "ion1", value = 0)),
                          column(3, numericInput("peion2", "ion2", value = 0))
                        ),
                        fluidRow(
                          column(3, verbatimTextOutput("pesn1")),
                          column(3, verbatimTextOutput("pesn2")),
                          column(3, verbatimTextOutput("pesum"))
                        )
               )
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for PEs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + H - phosphoethanolamina]+")),
        column(3, 
               strong("Negative [M-H]-:"),
               tags$li("[M - H - (sn1-H2O)]-"),
               tags$li("[sn1 - H]-"),
               tags$li("[M - H - (sn2-H2O)]-"),
               tags$li("[sn2 - H]-")
               
        )
      ), # close tab PEs
      
      ### Lyso-PC ----
      tabPanel(
        "Lyso-Phosphatidylcholines (Lyso-PCs)",
        h1("Lyso-Phosphatidylcholines (Lyso-PCs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("lpcC", "C", value = 16)),
                 column(2, numericInput("lpcdb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("lpcformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("lpcmzvals1")),
               fluidRow(verbatimTextOutput("lpcmzvals2"))
        ),
        column(1), 
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("lpcfragpos")),
               fluidRow(h4("ESI-"), verbatimTextOutput("lpcfragneg"))
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for Lyso-PCs:"),
        column(3,
               strong("Positive [M+H]+:"),
               tags$li("[M + H - H2O]+"),
               tags$li("[Phosphocholine]+")),
        column(3, 
               strong("Negative [M-H+HCOOH]-:"),
               tags$li("[M - H - CH3]-"),
               br(),
               strong("[M - H - CH3]-"),
               tags$li("[sn - H]-")
               
        )
      ), # close tab Lyso-PCs
      
      ### PC ----
      tabPanel(
        "Phosphatidylcholines (PCs)",
        h1("Phosphatidylcholines (PCs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("pcC", "C", value = 32)),
                 column(2, numericInput("pcdb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("pcformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("pcmzvals1")),
               fluidRow(verbatimTextOutput("pcmzvals2"))
        ),
        column(1), 
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("pcfragpos")),
               fluidRow(h4("ESI-"), verbatimTextOutput("pcfragneg"))
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for PCs:"),
        column(3,
               strong("Positive [M+Na]+:"),
               tags$li("[M + Na - C3H9N]+"),
               tags$li("[M + Na - phosphocholine]+")
        ),
        column(3, 
               strong("Negative [M-H+HCOOH]-:"),
               tags$li("[M - H - CH2]-")
        )
      ), # close tab PCs
      
      ### PI ----
      tabPanel(
        "Phosphatidylinositols (PIs)",
        h1("Phosphatidylinositols (PIs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("piC", "C", value = 32)),
                 column(2, numericInput("pidb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("piformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("pimzvals1")),
               fluidRow(verbatimTextOutput("pimzvals2"))
        ),
        column(1), 
        column(6, h3("MS2 (-)"),
               fluidRow(
                 column(3, numericInput("piion1", "ion1", value = 0)),
                 column(3, numericInput("piion2", "ion2", value = 0))
               ),
               fluidRow(
                 column(3, verbatimTextOutput("pisn1")),
                 column(3, verbatimTextOutput("pisn2")),
                 column(3, verbatimTextOutput("pisum"))
               ),
               hr(),
               fluidRow(
                 column(4, selectInput("piion1x", "sn1", choices = sn.list)),
                 column(4, selectInput("piion2x", "sn2", choices = sn.list))
               ),
               fluidRow(
                 column(4, verbatimTextOutput("pisn1x")),
                 column(4, verbatimTextOutput("pisn2x"))
               )
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for PIs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + H]+"),
               tags$li("[M + NH4 - phosphoinositol]+")),
        column(3, 
               strong("Negative [M-H]-:"),
               tags$li("[M - H - sn1]-"),
               tags$li("[M - H - (sn1-H2O)]-"),
               tags$li("[sn1 - H]-"),
               tags$li("[M - H - sn2]-"),
               tags$li("[M - H - (sn2-H2O)]-"),
               tags$li("[sn2 - H]-")
               
        )
      ) # close tab PIs
      
      
    ), # close Glycerophospholipids [GP]
    
    ## GL ----
    navbarMenu(
      "Glycerolipids (GL) - Glycerols",
      ### DAG ----
      tabPanel(
        "Diacylglycerols (DAGs)",
        h1("Diacylglycerols (DAGs)"),
        column(4, h3("Formula"),
               fluidRow(
                 column(2, numericInput("dagC", "C", value = 39)),
                 column(2, numericInput("dagdb", "db", value = 0))
               ),
               column(4, fluidRow(verbatimTextOutput("dagformula"))),
               fluidRow(),
               fluidRow(h3("m/z values"), verbatimTextOutput("dagmzvals1")),
               fluidRow("", verbatimTextOutput("dagmzvals2"))), 
        column(1), 
        column(6, h3("MS2 (+)"),
               fluidRow(verbatimTextOutput("dagfragpos")),
               fluidRow(
                 column(3, numericInput("dagion1", "ion1", value = 0)),
                 column(3, numericInput("dagion2", "ion2", value = 0))
               ),
               fluidRow(
                 column(3, verbatimTextOutput("dagsn1")),
                 column(3, verbatimTextOutput("dagsn2")),
                 column(3, verbatimTextOutput("dagsum"))
               )),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for DAGs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + H]+"),
               tags$li("[M + H - H2O]+"),
               tags$li("[M + H - sn1]+"))
      ), # close tab DAG
      
      ### TAG ----
      tabPanel(
        "Triacylglycerols (TAGs)",
        h1("Triacylglycerols (TAGs)"),
        column(2, h3("Formula"),
               fluidRow(
                 column(6, numericInput("tagC", "C", value = 54)),
                 column(6, numericInput("tagdb", "db", value = 0))
               ),
               fluidRow(verbatimTextOutput("tagformula")),
               fluidRow(h3("m/z values"), verbatimTextOutput("tagmzvals1")),
               fluidRow(verbatimTextOutput("tagmzvals2"))
        ),
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
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for TAGs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + H]+"),
               tags$li("[M + H - sn1]+"),
               tags$li("[M + H - sn2]+"),
               tags$li("[M + H - sn3]+")
        )
      ) # close tab TAG
    ), # close GL-glycerols
    navbarMenu(
      "Glycerolipids (GL) - Galactosyldiacylglycerols (GDG)",
      
      ### MGDG ----
      tabPanel(
        "Monogalactosyldiacylglycerols (MGDG)",
        h1("Monogalactosyldiacylglycerols (MGDG)"),
        column(3, h3("Formula"),
               fluidRow(
                 column(6, numericInput("mgdgC", "C", value = 36)),
                 column(6, numericInput("mgdgdb", "db", value = 0))
               ),
               fluidRow(verbatimTextOutput("mgdgformula")),
               fluidRow(h3("m/z values"), verbatimTextOutput("mgdgmzvals1")),
               fluidRow(verbatimTextOutput("mgdgmzvals2"))
        ),
        column(1),
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("mgdgfragpos")),
               fluidRow(column(2, strong("[M+Na]+")), column(3, numericInput("mgdgion1pos", "ion1", value = 0)),
                        column(3, numericInput("mgdgion2pos", "ion2", value = 0))),
               fluidRow(
                 column(2, ""),
                 column(3, verbatimTextOutput("mgdgsn1pos")),
                 column(3, verbatimTextOutput("mgdgsn2pos")),
                 column(3, verbatimTextOutput("mgdgsumpos"))
               ),
               fluidRow(h4("ESI-"), 
                        column(4, verbatimTextOutput("mgdgfragneg")),
                        column(2, numericInput("mgdgion1", "ion1", value = 0)),
                        column(6, verbatimTextOutput("mgdgsn1"))
               )
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for MGDGs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + NH4 - galactose]+ + 0.984"),
               tags$li("[M + H - galactose]+"),
               br(),
               strong("[M+Na]+:"),
               tags$li("[M + Na - sn1]+"),
               tags$li("[M + NA - sn2]+")),
        column(3, 
               strong("Negative [M-H]-:"),
               tags$li("[M - H]-"),
               tags$li("[M - H - (sn1 - H2O)]-"),
               tags$li("[sn1 - H]-")
        )
      ), # close MGDG
      
      ### DGDG ----
      tabPanel(
        "Digalactosyldiacylglycerol (DGDG)",
        h1("Digalactosyldiacylglycerol (DGDG)"),
        column(2, h3("Formula"),
               fluidRow(
                 column(6, numericInput("dgdgC", "C", value = 36)),
                 column(6, numericInput("dgdgdb", "db", value = 0))
               ),
               fluidRow(verbatimTextOutput("dgdgformula")),
               fluidRow(h3("m/z values"), verbatimTextOutput("dgdgmzvals1")),
               fluidRow(verbatimTextOutput("dgdgmzvals2"))
        ),
        column(1),
        column(6, h3("MS2"),
               fluidRow(h4("ESI+"), verbatimTextOutput("dgdgfragpos")),
               fluidRow(h4("ESI-"), 
                        column(6, verbatimTextOutput("dgdgfragneg")),
                        column(3, numericInput("dgdgion1", "ion1", value = 0)),
                        column(3, verbatimTextOutput("dgdgsn1"))
               )
        ),
        fluidRow(),
        hr(),
        h3("Commonly occuring product ions for DGDGs:"),
        column(3,
               strong("Positive [M+NH4]+:"),
               tags$li("[M + NH4 - 2(galactose)]+ + 0.984"),
               tags$li("[M + H - 2(galactose)]+")
        ),
        column(3, 
               strong("Negative [M-H]-:"),
               tags$li("[M - H]-"),
               tags$li("[M - H - sn1]-")
        )
      ) # close DGDG
    ), #  close Glycosylglycerols
    
    ## Carnitines ----
    tabPanel(
      "Carnitines",
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("carC", "C", value = 24)),
               column(2, numericInput("cardb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("carformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("carmzvals1")),
             fluidRow(verbatimTextOutput("carmzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("carfragpos")),
             fluidRow(h4("ESI-"), verbatimTextOutput("carfragneg"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for Carnitines:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - H2O]+"),
             tags$li("[Phosphocholine]+")),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - H - CH3]-"),
             br(),
             strong("[M - H - CH3]-"),
             tags$li("[sn - H]-")
             
      )
    ) # close tab Carnitines
    
  )# close ui
  
  # SERVER ---------------------------------------------------------------------
  server <- function(input, output) {
    
    output$posfrag1add <- renderPrint({
      idx <- c(which(abs((input$posprec - input$posfrag1) - mzdif.pos$dif) < 0.01),
               unlist(matchWithPpm(input$posfrag1, mzdif.pos$dif, ppm = 10)))
      mzdif.pos$add[idx]
    })
    
    output$posfrag2add <- renderPrint({
      idx <- c(which(abs((input$posprec - input$posfrag2) - mzdif.pos$dif) < 0.01),
               unlist(matchWithPpm(input$posfrag2, mzdif.pos$dif, ppm = 10)))
      mzdif.pos$add[idx]
    })
    
    output$posfrag3add <- renderPrint({
      idx <- c(which(abs((input$posprec - input$posfrag3) - mzdif.pos$dif) < 0.01),
               unlist(matchWithPpm(input$posfrag3, mzdif.pos$dif, ppm = 10)))
      mzdif.pos$add[idx]
    })
    
    output$posfrag4add <- renderPrint({
      idx <- c(which(abs((input$posprec - input$posfrag4) - mzdif.pos$dif) < 0.01),
               unlist(matchWithPpm(input$posfrag4, mzdif.pos$dif, ppm = 10)))
      mzdif.pos$add[idx]
    })
    
    output$negfrag1add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag1) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag1, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    output$negfrag2add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag2) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag2, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    output$negfrag3add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag3) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag3, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    output$negfrag4add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag4) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag4, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    output$negfrag5add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag5) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag5, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    output$negfrag6add <- renderPrint({
      idx <- c(which(abs((input$negprec - input$negfrag6) - mzdif.neg$dif) < 0.01),
               unlist(matchWithPpm(input$negfrag6, mzdif.neg$dif, ppm = 10)))
      mzdif.neg$add[idx]
    })
    
    ## PA ----
    pafml <- reactive({
      paste0("C", input$paC + 3, "H", 
             input$paC*2 - (2 + 2*input$padb) + 7, "O8P")
    })
    
    pamass <- reactive({
      MonoisotopicMass(formula = ListFormula(pafml()))
    })
    
    output$paformula <- renderPrint({pafml()})
    
    output$pamzvals1 <- renderPrint({
      mass2mz(pamass(), adduct = c("[M+NH4]+", "[M-H]-"))
    })
    
    output$pamzvals2 <- renderPrint({
      mass2mz(
        pamass(), 
        adduct = c("[2M+H]+", "[2M+NH4]+", "[2M-H]-"))
    })
    
    output$pafragpos <- renderPrint({
      tmp <- mass2mz(pamass(),  adduct = c("[M+H]+"))
      tmp2 <- colnames(tmp)
      tmp <- c(as.numeric(unlist(mass2mz(pamass(), "[M+H]+"))) - 
                 MonoisotopicMass(formula = ListFormula("H3PO4")), tmp)
      names(tmp) <- c("[M+H-PA]+", tmp2)
      tmp
    })
    
    output$pasn1 <- renderPrint({
      idx <- unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion1, 
        sn$mass, ppm = 10))
      HTML(paste(
        sn$sn[idx], 
        round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - sn$mass[idx] + MonoisotopicMass(formula = ListFormula("H2O")), 4), 
        round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - sn$mass[idx], 4), 
        round(mass2mz(sn$mass[idx], "[M-H]-"), 4),
        sep = '<br/>'
      ))
    })
    
    output$pasn2 <- renderPrint({
      idx <- unlist(matchWithPpm(
        unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - input$paion2, 
        sn$mass, ppm = 10))
      HTML(paste(
        sn$sn[idx], 
        round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - sn$mass[idx] + MonoisotopicMass(formula = ListFormula("H2O")), 4), 
        round(unlist(mass2mz(pamass(), adduct = c("[M-H]-"))) - sn$mass[idx], 4), 
        round(mass2mz(sn$mass[idx], "[M-H]-"), 4),
        sep = '<br/>'
      ))
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
    
    output$mpamzvals1 <- renderPrint({
      mass2mz(mpamass(), adduct = c("[M+NH4]+", "[M-H]-"))
    })
    
    output$mpamzvals2 <- renderPrint({
      mass2mz(mpamass(), adduct = c("[2M+H]+", "[2M+NH4]+", "[2M-H]-"))
    })
    
    output$mpafragpos <- renderPrint({
      tmp <- mass2mz(mpamass(), "[M+H]+")
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
    
    output$mpasum <- renderPrint({
      paste0(
        sn$C[unlist(matchWithPpm(input$mpaion1 + 1.007276, sn$mass, ppm = 10))] +
          sn$C[unlist(matchWithPpm(input$mpaion2 + 1.007276, sn$mass, ppm = 10))],
        ":",
        sn$db[unlist(matchWithPpm(input$mpaion1 + 1.007276, sn$mass, ppm = 10))] +
          sn$db[unlist(matchWithPpm(input$mpaion2 + 1.007276, sn$mass, ppm = 10))]
      )
    })
    
    ## PG ----
    pgfml <- reactive({
      paste0("C", input$pgC + 6, "H", 
             input$pgC*2 - (2 + 2*input$pgdb) + 13, "O10P")
    })
    
    pgmass <- reactive({
      MonoisotopicMass(formula = ListFormula(pgfml()))
    })
    
    output$pgformula <- renderPrint({pgfml()})
    
    output$pgmzvals1 <- renderPrint({
      mass2mz(pgmass(), adduct = c("[M+NH4]+", "[M-H]-"))
    })
    
    #output$pgmzvals2 <- renderPrint({
    #  tmp <- unlist(mass2mz(
    #    pgmass(), 
    #    adduct = c("[M+NH4]+", "[2M+NH4]+"#, ""[2M-H]-"
    #    )))
    #  tmp2 <- colnames(tmp)
    #  tmp <- c(as.numeric(unlist(mass2mz(pgmass(), "[M+NH4]+"))) - 
    #             MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984,
    #           as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
    #             MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
    #  names(tmp) <- c("[M+H-pg]+", "[M+C2H8N]+", tmp2)
    #  tmp <- tmp[1, 3, 2, 4]
    #  tmp
    #})
    
    #output$pgfragpos <- renderPrint({
    #  round(as.numeric(mass2mz(pgmass(), "[M+H]+")) - 
    #          MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 5)
    #})
    
    output$pgsn1 <- renderPrint({
      idx1 <- unlist(matchWithPpm(
        mass2mz(pgmass(), "[M-H]-") - input$pgion1, 
        c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
      idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion1, sn$mass, ppm = 10))
      idx <- c(idx1, idx2)
      HTML(paste(sn$sn[idx], 
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - sn$mass[idx]),
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (
                   (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) +
                     MonoisotopicMass(formula = ListFormula("C3H8O3")))),
                 sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
                 sep = '<br/>'))
    })
    
    output$pgsn2 <- renderPrint({
      idx1 <- unlist(matchWithPpm(
        mass2mz(pgmass(), "[M-H]-") - input$pgion2, 
        c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
      idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion2, sn$mass, ppm = 10))
      idx <- c(idx1, idx2)
      HTML(paste(sn$sn[idx], 
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - sn$mass[idx]),
                 sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (
                   (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) +
                     MonoisotopicMass(formula = ListFormula("C3H8O3")))),
                 sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
                 sep = '<br/>'))
    })
    
    output$pgsum <- renderPrint({
      idx1 <- unlist(matchWithPpm(
        mass2mz(pgmass(), "[M-H]-") - input$pgion1, 
        c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
      idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion1, sn$mass, ppm = 10))
      idxa <- c(idx1, idx2)
      idx1 <- unlist(matchWithPpm(
        mass2mz(pgmass(), "[M-H]-") - input$pgion2, 
        c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
      idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion2, sn$mass, ppm = 10))
      idxb <- c(idx1, idx2)
      paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
    })
    
    ## PE ----
    pefml <- reactive({
      paste0("C", input$peC + 5, "H", 
             input$peC*2 - (2 + 2*input$pedb) + 12, "NO8P")
    })
    
    pemass <- reactive({
      MonoisotopicMass(formula = ListFormula(pefml()))
    })
    
    output$peformula <- renderPrint({pefml()})
    
    output$pemzvals1 <- renderPrint({
      mass2mz(pemass(), adduct = c("[M+H]+", "[M-H]-"))
    })
    
    #output$pemzvals2 <- renderPrint({
    #  tmp <- unlist(mass2mz(
    #    pemass(), 
    #    adduct = c("[M+NH4]+", "[2M+NH4]+"#, ""[2M-H]-"
    #    )))
    #  tmp2 <- colnames(tmp)
    #  tmp <- c(as.numeric(unlist(mass2mz(pemass(), "[M+NH4]+"))) - 
    #             MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984,
    #           as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
    #             MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
    #  names(tmp) <- c("[M+H-pe]+", "[M+C2H8N]+", tmp2)
    #  tmp <- tmp[1, 3, 2, 4]
    #  tmp
    #})
    
    output$pefragpos <- renderPrint({
      round(as.numeric(mass2mz(pemass(), "[M+H]+")) - 
              MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 5)
    })
    
    output$pesn1 <- renderPrint({
      paste0(
        sn$sn[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion1, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))], " (",
        round(mass2mz(sn$mass[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion1, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))], "[M-H]-"), 5), ")"
      )
    })
    
    output$pesn2 <- renderPrint({
      paste0(
        sn$sn[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion2, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))], " (",
        round(mass2mz(sn$mass[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion2, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))], "[M-H]-"), 5), ")"
      )
    })
    
    output$pesum <- renderPrint({
      paste0(
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion1, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))] +
          sn$C[unlist(matchWithPpm(
            unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion2, 
            sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))],
        ":",
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion1, 
          sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))] +
          sn$db[unlist(matchWithPpm(
            unlist(mass2mz(pemass(), adduct = c("[M-H]-"))) - input$peion2, 
            sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))]
      )
    })
    
    ## Lyso-PC ----
    lpcfml <- reactive({
      paste0("C", input$lpcC + 8, "H", 
             input$lpcC*2 - (2 + 2*input$lpcdb) + 20, "NO7P")
    })
    
    lpcmass <- reactive({
      MonoisotopicMass(formula = ListFormula(lpcfml()))
    })
    
    output$lpcformula <- renderPrint({lpcfml()})
    
    output$lpcmzvals1 <- renderPrint({
      mass2mz(lpcmass(), adduct = c("[M+H]+", "[M+CHO2]-"))
    })
    
    #output$lpcmzvals2 <- renderPrint({
    #  tmp <- unlist(mass2mz(
    #    lpcmass(), 
    #    adduct = c("[M+NH4]+", "[2M+NH4]+"#, ""[2M-H]-"
    #    )))
    #  tmp2 <- colnames(tmp)
    #  tmp <- c(as.numeric(unlist(mass2mz(lpcmass(), "[M+NH4]+"))) - 
    #             MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984,
    #           as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
    #             MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
    #  names(tmp) <- c("[M+H-lpc]+", "[M+C2H8N]+", tmp2)
    #  tmp <- tmp[1, 3, 2, 4]
    #  tmp
    #})
    
    output$lpcfragpos <- renderPrint({
      c(round(MonoisotopicMass(formula = ListFormula("C5H15NO4P")), 5),
        round(mass2mz(lpcmass(), "[M+H-H2O]+"), 5))
    })
    
    output$lpcfragneg <- renderPrint({
      c(round(mass2mz(lpcmass(), "[M-H]-") - MonoisotopicMass(formula = ListFormula("CH2")), 5),
        mass2mz(sn$mass[sn$sn == paste0(input$lpcC, ":", input$lpcdb)], "[M-H]-"))
    })
    
    ## PC ----
    pcfml <- reactive({
      paste0("C", input$pcC + 8, "H", 
             input$pcC*2 - (2 + 2*input$pcdb) + 18, "NO8P")
    })
    
    pcmass <- reactive({
      MonoisotopicMass(formula = ListFormula(pcfml()))
    })
    
    output$pcformula <- renderPrint({pcfml()})
    
    output$pcmzvals1 <- renderPrint({
      mass2mz(pcmass(), adduct = c("[M+H]+", "[M+CHO2]-"))
    })
    
    output$pcmzvals2 <- renderPrint({
      mass2mz(pcmass(), adduct = c("[M+Na]+", "[2M+H]+", "[2M+CHO2]-"))
    })
    
    output$pcfragpos <- renderPrint({
      c(as.numeric(mass2mz(pcmass(), "[M+Na]+")) - MonoisotopicMass(formula = ListFormula("C3H9N")),
        as.numeric(mass2mz(pcmass(), "[M+Na]+")) - MonoisotopicMass(formula = ListFormula("C5H14NO4P")))
    })
    
    output$pcfragneg <- renderPrint({
      as.numeric(mass2mz(pcmass(), "[M-H]-")) - MonoisotopicMass(formula = ListFormula("CH2"))
    })
    
    ## PI ----
    pifml <- reactive({
      paste0("C", input$piC + 9, "H", 
             input$piC*2 - (2 + 2*input$pidb) + 17, "O13P")
    })
    
    pimass <- reactive({
      MonoisotopicMass(formula = ListFormula(pifml()))
    })
    
    output$piformula <- renderPrint({pifml()})
    
    output$pimzvals1 <- renderPrint({
      mass2mz(pimass(), adduct = c("[M+NH4]+", "[M-H]-"))
    })
    
    output$pimzvals2 <- renderPrint({
      tmp <- unlist(mass2mz(
        pimass(), 
        adduct = c("[M+H]+", "[2M+NH4]+"#, ""[2M-H]-"
        )))
      tmp2 <- colnames(tmp)
      tmp <- c(as.numeric(unlist(mass2mz(pimass(), "[M+NH4]+"))) - 
                 MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984,
               as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
                 MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
      names(tmp) <- c("[M+H-PI]+", "[M+C2H8N]+", tmp2)
      tmp <- tmp[1, 3, 2, 4]
      tmp
    })
    
    output$pisn1 <- renderPrint({
      paste0(
        sn$sn[unlist(matchWithPpm(
          unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion1, 
          sn$mass, ppm = 10))], " (",
        round(unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - 
                sn$mass[unlist(matchWithPpm(
                  unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion1, 
                  sn$mass, ppm = 10))] + 
                MonoisotopicMass(formula = ListFormula("H2O")), 4), ")"
      )
    })
    
    output$pisn2 <- renderPrint({
      paste0(
        sn$sn[unlist(matchWithPpm(
          unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion2, 
          sn$mass, ppm = 10))], " (",
        round(unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - 
                sn$mass[unlist(matchWithPpm(
                  unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion2, 
                  sn$mass, ppm = 10))] + 
                MonoisotopicMass(formula = ListFormula("H2O")), 4), ")"
      )
    })
    
    output$pisum <- renderPrint({
      paste0(
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion1, 
          sn$mass, ppm = 10))] +
          sn$C[unlist(matchWithPpm(
            unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion2, 
            sn$mass, ppm = 10))],
        ":",
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion1, 
          sn$mass, ppm = 10))] +
          sn$db[unlist(matchWithPpm(
            unlist(mass2mz(pimass(), adduct = c("[M-H]-"))) - input$piion2, 
            sn$mass, ppm = 10))]
      )
    })
    
    output$pisn1x <- reactive({
      round(unlist(mass2mz(sn$mass[sn$sn == input$piion1x], "[M-H]-")), 5)
    })
    
    output$pisn2x <- reactive({
      round(unlist(mass2mz(sn$mass[sn$sn == input$piion2x], "[M-H]-")), 5)
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
    
    output$dagmzvals1 <- renderPrint({
      unlist(mass2mz(
        dagmass(), 
        adduct = c("[M+NH4]+", "[M+CHO2]-")))
    })
    
    output$dagmzvals2 <- renderPrint({
      tmp <- unlist(mass2mz(
        dagmass(), 
        adduct = c("[M+Na]+", "[2M+NH4]+", "[2M+Na]+", "[2M+CHO2]-")))
      tmp2 <- colnames(tmp)
      tmp <- c(tmp, as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
                 MonoisotopicMass(formula = ListFormula("C2H7N")))
      names(tmp) <- c(tmp2, "[M+C2H8N]+")
      tmp <- tmp[c(1:3, 7, 4:6)]
      tmp
    })
    
    output$dagfragpos <- renderPrint({
      mass2mz(
        dagmass(), 
        adduct = c("[M+H-H2O]+", "[M+H]+"))
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
    
    output$tagmzvals1 <- renderPrint({
      mass2mz(tagmass(), "[M+NH4]+")
    })
    
    output$tagmzvals2 <- renderPrint({
      mass2mz(tagmass(), "[M+H]+")
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
    
    ## MGDG -----
    
    mgdgfml <- reactive({
      paste0("C", input$mgdgC + 9, "H", 
             input$mgdgC*2 - (2 + 2*input$mgdgdb) + 16, "O10")
    })
    
    mgdgmass <- reactive({
      MonoisotopicMass(formula = ListFormula(mgdgfml()))
    })
    
    output$mgdgformula <- renderPrint({mgdgfml()})
    
    output$mgdgmzvals1 <- renderPrint({
      unlist(mass2mz(mgdgmass(), adduct = c("[M+NH4]+", "[M+CHO2]-")))
    })
    
    output$mgdgmzvals2 <- renderPrint({
      tmp <- unlist(mass2mz(
        mgdgmass(), 
        adduct = c("[M+H]+", "[M+Na]+", "[2M+NH4]+", "[2M+Na]+")))
      tmp2 <- colnames(tmp)
      tmp <- c(as.numeric(unlist(mass2mz(mgdgmass(), "[M+H]+"))) + 
                 MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
      names(tmp) <- c("[M+C2H8N]+", tmp2)
      tmp <- tmp[c(2:3, 1, 4:5)]
      tmp
    })
    
    output$mgdgfragpos <- renderPrint({
      paste(
        round(as.numeric(mass2mz(mgdgmass(), "[M+NH4]+")) - 
                MonoisotopicMass(formula = ListFormula("C6H12O6")) + 0.984, 5), "-",
        round(as.numeric(mass2mz(mgdgmass(), "[M+H]+")) - 
                MonoisotopicMass(formula = ListFormula("C6H12O6")), 5))
    })
    
    output$mgdgsn1pos <- renderPrint({
      sn$sn[unlist(matchWithPpm(
        unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion1pos, 
        sn$mass, ppm = 10))]
    })
    
    output$mgdgsn2pos <- renderPrint({
      sn$sn[unlist(matchWithPpm(
        unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion2pos, 
        sn$mass, ppm = 10))]
    })
    
    output$mgdgsumpos <- renderPrint({
      paste0(
        sn$C[unlist(matchWithPpm(
          unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion1pos, 
          sn$mass, ppm = 10))] +
          sn$C[unlist(matchWithPpm(
            unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion2pos, 
            sn$mass, ppm = 10))],
        ":",
        sn$db[unlist(matchWithPpm(
          unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion1pos, 
          sn$mass, ppm = 10))] +
          sn$db[unlist(matchWithPpm(
            unlist(mass2mz(mgdgmass(), adduct = c("[M+Na]+"))) - input$mgdgion2pos, 
            sn$mass, ppm = 10))]
      )
    })
    
    output$mgdgfragneg <- renderPrint({
      round(mass2mz(mgdgmass(), "[M-H]-"), 5)
    })
    
    output$mgdgsn1 <- renderPrint({
      idx <- unlist(matchWithPpm(
        unlist(mass2mz(mgdgmass(), adduct = c("[M-H]-"))) - input$mgdgion1, 
        sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), ppm = 10))
      paste0(sn$sn[idx], " (", round(mass2mz(sn$mass[idx], "[M-H]-"), 5), ") - ",
             input$mgdgC - sn$C[idx], ":", input$mgdgdb - sn$db[idx],
             " (",
             round(mass2mz(sn$mass[sn$sn == paste0(input$mgdgC - sn$C[idx], ":", 
                                                   input$mgdgdb - sn$db[idx])], "[M-H]-"), 5)
             , ")")
    })
    
    ## DGDG -----
    
    dgdgfml <- reactive({
      paste0("C", input$dgdgC + 15, "H", 
             input$dgdgC*2 - (2 + 2*input$dgdgdb) + 26, "O15")
    })
    
    dgdgmass <- reactive({
      MonoisotopicMass(formula = ListFormula(dgdgfml()))
    })
    
    output$dgdgformula <- renderPrint({dgdgfml()})
    
    output$dgdgmzvals1 <- renderPrint({
      unlist(mass2mz(dgdgmass(), adduct = c("[M+NH4]+", "[M+CHO2]-")))
    })
    
    output$dgdgmzvals2 <- renderPrint({
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+", "[M-H]-")))
    })
    
    output$dgdgfragpos <- renderPrint({
      paste(
        round(as.numeric(mass2mz(dgdgmass(), "[M+NH4]+")) - 
                MonoisotopicMass(formula = ListFormula("C12H22O11")) + 0.984, 5), "-",
        round(as.numeric(mass2mz(dgdgmass(), "[M+H]+")) - 
                MonoisotopicMass(formula = ListFormula("C12H22O11")), 5))
    })
    
    output$dgdgfragneg <- renderPrint({
      round(mass2mz(dgdgmass(), "[M-H]-"), 5)
    })
    
    output$dgdgsn1 <- renderPrint({
      idx <- unlist(matchWithPpm(
        unlist(mass2mz(dgdgmass(), adduct = c("[M-H]-"))) - input$dgdgion1, 
        sn$mass, ppm = 10))
      paste0(sn$sn[idx], "-",
             input$dgdgC - sn$C[idx], ":", input$dgdgdb - sn$db[idx])
    })
    
    ## Carnitines ----
    carfml <- reactive({
      paste0("C", input$carC + 7, "H", 
             input$carC*2 - (2 + 2*input$cardb) + 15, "NO4")
    })
    
    carmass <- reactive({
      MonoisotopicMass(formula = ListFormula(carfml()))
    })
    
    output$carformula <- renderPrint({carfml()})
    
    output$carmzvals1 <- renderPrint({
      mass2mz(carmass(), adduct = c("[M+H]+", "[M+CHO2]-"))
    })
    
    output$carfragpos <- renderPrint({
      round(as.numeric(mass2mz(carmass(), "[M+H]+")) - MonoisotopicMass(formula = ListFormula("C3H9N")), 5)
    })
    
    output$carfragneg <- renderPrint({
      round(as.numeric(mass2mz(carmass(), "[M-H]-")) - MonoisotopicMass(formula = ListFormula("C7H13NO2")), 5)
    })
    
  } # close server
  
  shinyApp(ui = ui, server = server)
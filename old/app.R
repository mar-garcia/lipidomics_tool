options(repos = BiocManager::repositories())
options("repos")

library(shiny)
library(MetaboCoreUtils)
library(OrgMassSpecR)
library(DT)
library(Spectra)

.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}
label_fun <- function(x) {
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 7)
  mzs[ints < 5] <- ""
  mzs
}

sn <- data.frame(
  "C" = rep(seq(14, 24), each = 4),
  "db" = rep(seq(0, 3), 11)
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
  c(MonoisotopicMass(formula = ListFormula("NH3")), 
    "loss NH3 -> PA / mPA / dmPA / PI / DAG / TAG / MGDG / DGDG"),
  c(MonoisotopicMass(formula = ListFormula("H2O")), 
    "loss H2O -> CER / Lyso PA / Lyso PE / Lyso PC / Lyso PS"),
  c(MonoisotopicMass(formula = ListFormula("NH3H2O")), 
    "loss NH3 & H2O -> DAG"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4NH3")), 
    "loss NH3 & phosphate -> PA / dmPA"),
  c(MonoisotopicMass(formula = ListFormula("NH3H3PO4CH2")), 
    "loss NH3 & mPA -> mPA"),
  c(MonoisotopicMass(formula = ListFormula("NH3H3PO4C2H4")), 
    "loss NH3 & dmPA -> dmPA"),
  c(MonoisotopicMass(formula = ListFormula("C3H9N")), 
    "loss trimethylamine -> PC / Carnitine"),
  c(MonoisotopicMass(formula = ListFormula("C5H14NO4P")), 
    "loss C5H14NO4P -> PC"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) - 0.984, 
    "loss phosphoinositol (NH4 ion) -> PI"),
  c(MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 
    "loss phosphoethanolamine -> Lyso PE / PE"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> MGDG [M+Na]+")),
  cbind(sn$mass + MonoisotopicMass(formula = ListFormula("NH3")), 
        paste("loss NH3 &", sn$sn, "-> DAG / TAG")),
  cbind(sn$mass, paste("loss ", sn$sn, "-> PC")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), 
        paste0("loss '", sn$sn, "-H2O' -> PC")),
  
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
  cbind(MonoisotopicMass(formula = ListFormula("C5H15NO4P")), 
        "[Phosphocholine]+ -> Lyso PC")
))
colnames(mzdif.pos) <- c("dif", "add")
mzdif.pos$dif <- as.numeric(mzdif.pos$dif)

mzdif.neg <- data.frame(rbind(
  c(MonoisotopicMass(formula = ListFormula("HCOOH")), 
    "loss HCOOH -> CER / MGDG / DGDG"),
  c(MonoisotopicMass(formula = ListFormula("CO2")), 
    "loss of CO2 (HCOOH ion) -> FFA"),
  c(MonoisotopicMass(formula = ListFormula("C3H5NO2")), 
    "loss C3H5NO2 -> LysoPS"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> Lyso PA / PA / PG / PI")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), 
        paste0("loss '", sn$sn, "-H2O' -> PA / mPA / dmPA / PG / PE")),
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
  c(MonoisotopicMass(formula = ListFormula("HCOOHC7H13NO2")), 
    "loss HCOOH & carnitine -> Carnitine"),
  cbind(mass2mz(sn$mass, "[M-H]-"), 
        paste0("[", sn$sn, 
               "-H]- -> PA / mPA / dmPA / PG / Lyso PG / PE / Lyso PC / PI / MGDG"))
))
colnames(mzdif.neg) <- c("dif", "add")
mzdif.neg$dif <- as.numeric(mzdif.neg$dif)


load("MS2_spectra.RData")
sps_ms2 <- c(sps_ms2_POS_cmps, sps_ms2_NEG_cmps)
sps_ms2_nl <- sps_ms2
sps_ms2_nl <- applyProcessing(sps_ms2_nl)
mz(sps_ms2_nl@backend) <- mz(sps_ms2_nl) - precursorMz(sps_ms2_nl)

sps_ms2_std <- c(sps_ms2_POS_stds, sps_ms2_NEG_stds)
sps_ms2_std_nl <- sps_ms2_std
sps_ms2_std_nl <- applyProcessing(sps_ms2_std_nl)
mz(sps_ms2_std_nl@backend) <- mz(sps_ms2_std_nl) - precursorMz(sps_ms2_std_nl)

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
               column(2, numericInput("paC", "C", value = 36)),
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
             br(),
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
               column(2, numericInput("mpaC", "C", value = 36)),
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
                        column(3, htmlOutput("mpasn1")),
                        column(3, htmlOutput("mpasn2")),
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
             tags$li("[M - H - (sn1-H2O)]-"),
             tags$li("[sn1 - H]-"),
             br(),
             tags$li("[M - H - (sn2-H2O)]-"),
             tags$li("[sn2 - H]-"))
    ), # close mPA
    
    ### PA dimethylated ----
    tabPanel(
      "Dimethylated phosphatidic acids (dmPAs)",
      h1("Methylated phosphatidic acids (dmPAs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("dmpaC", "C", value = 36)),
               column(2, numericInput("dmpadb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("dmpaformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("dmpamzvals1")),
             fluidRow(verbatimTextOutput("dmpamzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("dmpafragpos")),
             fluidRow(h4("ESI-"), 
                      fluidRow(
                        column(3, numericInput("dmpaion1", "ion1", value = 0)),
                        column(3, numericInput("dmpaion2", "ion2", value = 0))
                      ),
                      fluidRow(
                        column(3, htmlOutput("dmpasn1")),
                        column(3, htmlOutput("dmpasn2")),
                        column(3, verbatimTextOutput("dmpasum"))
                      ))),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for dmpas:"),
      column(3,
             strong("Positive [M+NH4]+:"),
             tags$li("[M + H]+"),
             tags$li("[M + H - phosphate]+"),
             tags$li("[M + H - dimethyl-phosphate]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[M - H - (sn1-H2O)]-"),
             tags$li("[sn1 - H]-"),
             tags$li("[sn2 - H]-"))
    ), # close dmPA
    
    ### PG ----
    tabPanel(
      "Phosphatidylglycerols (PGs)",
      h1("Phosphatidylglycerols (PGs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("pgC", "C", value = 36)),
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
             strong("Positive [M+NH4]+:"),
             tags$li("[M + NH4 - 141.04]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[sn1 - H]-"),
             tags$li("[M - H - sn1]-"),
             tags$li("[M - H - sn1 - H2O]-"),
             tags$li("[M - H - sn1 - H2O - glycerol]-"),
             br(),
             tags$li("[sn2 - H]-"),
             tags$li("[M - H - sn2]-"),
             tags$li("[M - H - sn2 - H2O]-"),
             tags$li("[M - H - sn2 - H2O - glycerol]-")
             
      )
    ), # close tab PGs
    
    ### Lyso-PE ----
    tabPanel(
      "Lyso-Phosphatidylethanolamines (Lyso-PEs)",
      h1("Lyso-Phosphatidylethanolamines (Lyso-PEs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("lpeC", "C", value = 18)),
               column(2, numericInput("lpedb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("lpeformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("lpemzvals"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("lpefragpos")),
             fluidRow(h4("ESI-"), 
                      fluidRow(numericInput("lpeion1", "ion", value = 0)),
                      fluidRow(htmlOutput("lpesn")))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for Lyso-PEs:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - H2O]+"),
             tags$li("[M + H - phosphoethanolamina]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[sn - H]-")
      )
    ), # close tab Lyso-PE
    
    ### PE ----
    tabPanel(
      "Phosphatidylethanolamines (PEs)",
      h1("Phosphatidylethanolamines (PEs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("peC", "C", value = 36)),
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
                        column(3, htmlOutput("pesn1")),
                        column(3, htmlOutput("pesn2")),
                        column(3, verbatimTextOutput("pesum"))
                      )
             )
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for PEs:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - phosphoethanolamina]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[sn1 - H]-"),
             tags$li("[sn2 - H]-"),
             tags$li("[M - H - sn2]-")
             
      )
    ), # close tab PEs
    
    ### Lyso-PC ----
    tabPanel(
      "Lyso-Phosphatidylcholines (Lyso-PCs)",
      h1("Lyso-Phosphatidylcholines (Lyso-PCs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("lpcC", "C", value = 18)),
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
             tags$li("[M - CH3]-"),
             br(),
             strong("[M - CH3]-"),
             tags$li("[sn - H]-")
             
      )
    ), # close tab Lyso-PCs
    
    ### PC ----
    tabPanel(
      "Phosphatidylcholines (PCs)",
      h1("Phosphatidylcholines (PCs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(3, numericInput("pcC", "C", value = 36)),
               column(3, numericInput("pcdb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("pcformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("pcmzvals1")),
             fluidRow(verbatimTextOutput("pcmzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), 
                      fluidRow(strong("[M+Na]+:"), 
                               verbatimTextOutput("pcfragposna"))),
             br(),
             fluidRow(
               strong("[M+H]+:"),
               fluidRow(
                 verbatimTextOutput("pcfragposh")
               ),
               fluidRow(
                 column(3, numericInput("pcion1", "ion1", value = 0)),
                 column(3, numericInput("pcion2", "ion2", value = 0))
               ),
             ),
             fluidRow(
               column(3, htmlOutput("pcsn1")),
               column(3, htmlOutput("pcsn2")),
               column(3, verbatimTextOutput("pcsum"))
             ),
             fluidRow(h4("ESI-"), verbatimTextOutput("pcfragneg"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for PCs:"),
      column(3,
             strong("Positive [M+Na]+:"),
             tags$li("[M + Na - trimethylamine]+"),
             tags$li("[M + Na - phosphocholine]+"),
             br(),
             strong("[M+H]+:"),
             tags$li("< C34 [M + H - 157.0504]+ ; > C36 [M + H - 183.0684]+"),
             tags$li("[M + H - sn1]+"),
             tags$li("[M + H - (sn1 - H2O)]+"),
             tags$li("[M + H - sn2]+"),
             tags$li("[M + H - (sn2 - H2O)]+")
      ),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - CH3]-")
      )
    ), # close tab PCs
    
    ### Lyso-PS ----
    tabPanel(
      "Lyso-Phosphatidylserines (Lyso-PSs)",
      h1("Lyso-Phosphatidylserines (Lyso-PSs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("lpsC", "C", value = 18)),
               column(2, numericInput("lpsdb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("lpsformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("lpsmzvals1")),
             fluidRow(verbatimTextOutput("lpsmzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("lpsfragpos")),
             fluidRow(h4("ESI-"), verbatimTextOutput("lpsfragneg"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for Lyso-PCs:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - H2O]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[M - H - C3H5NO2]-")
             
      )
    ), # close tab Lyso-PSs
    
    
    
    ### PS ----
    tabPanel(
      "Phosphatidylserines (PSs)",
      h1("Phosphatidylserines (PSs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("psC", "C", value = 40)),
               column(2, numericInput("psdb", "db", value = 2))
             ),
             column(4, fluidRow(verbatimTextOutput("psformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("psmzvals"))),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("psfragpos")),
             fluidRow(h4("ESI-"),
                      fluidRow(verbatimTextOutput("psfragneg")), 
                      fluidRow(
                        column(3, numericInput("psion1", "ion1", value = 0))),
                      fluidRow(
                        column(3, htmlOutput("pssn1")),
                        column(3, htmlOutput("pssn2"))
                      )
             )),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for PSs:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - phosphoserine]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[M - H - serine]-"),
             tags$li("[M - H - serine - sn1]-")
      )
    ), # close tab PSs
    
    ### Lyso-PI ----
    tabPanel(
      "Lyso-Phosphatidylserines (Lyso-PIs)",
      h1("Lyso-Phosphatidylserines (Lyso-PIs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("lpiC", "C", value = 18)),
               column(2, numericInput("lpidb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("lpiformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("lpimzvals1")),
             fluidRow(verbatimTextOutput("lpimzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("lpifragpos")),
             fluidRow(h4("ESI-"), verbatimTextOutput("lpifragneg"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for Lyso-PCs:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[]-")
             
      )
    ), # close tab Lyso-PIs
    
    ### PI ----
    tabPanel(
      "Phosphatidylinositols (PIs)",
      h1("Phosphatidylinositols (PIs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("piC", "C", value = 36)),
               column(2, numericInput("pidb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("piformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("pimzvals1")),
             fluidRow(verbatimTextOutput("pimzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("pifragpos")),
             fluidRow(h4("ESI-"), 
                      fluidRow(
                        column(3, numericInput("piion1", "ion1", value = 0)),
                        column(3, numericInput("piion2", "ion2", value = 0))
                      ),
                      fluidRow(
                        column(3, htmlOutput("pisn1")),
                        column(3, htmlOutput("pisn2")),
                        column(3, verbatimTextOutput("pisum"))
                      )
             )),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for PIs:"),
      column(3,
             strong("Positive [M+NH4]+:"),
             tags$li("[M + H]+"),
             tags$li("[M + H - phosphoinositol]+")),
      column(3, 
             strong("Negative [M-H]-:"),
             tags$li("[sn1 - H]-"),
             tags$li("[M - H - sn1]-"),
             tags$li("[M - H - sn1 - H2O]-"),
             tags$li("[M - H - sn1 - H2O - inositol]-"),
             br(),
             tags$li("[sn2 - H]-"),
             tags$li("[M - H - sn2]-"),
             tags$li("[M - H - (sn2-H2O)]-"),
             tags$li("[M - H - (sn2-H2O) - inositol]-")
             
      )
    ) # close tab PIs
    
    
  ), # close Glycerophospholipids [GP]
  
  ## GL ----
  navbarMenu(
    "Glycerolipids (GL)",
    ### DAG ----
    tabPanel(
      "Diacylglycerols (DAGs)",
      h1("Diacylglycerols (DAGs)"),
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("dagC", "C", value = 36)),
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
    ), # close tab TAG
    
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
             fluidRow(h4("ESI+"), 
                      column(2, strong("[M+NH4]+")),
                      column(6, verbatimTextOutput("mgdgfragpos"))),
             fluidRow(column(2, strong("[M+Na]+")), 
                      column(3, numericInput("mgdgion1pos", "ion1", value = 0)),
                      column(3, numericInput("mgdgion2pos", "ion2", value = 0))),
             fluidRow(
               column(2, ""),
               column(3, verbatimTextOutput("mgdgsn1pos")),
               column(3, verbatimTextOutput("mgdgsn2pos")),
               column(3, verbatimTextOutput("mgdgsumpos"))
             ),
             fluidRow(h4("ESI-"), verbatimTextOutput("mgdgfragneg")),
             fluidRow(
               column(2, numericInput("mgdgion1", "ion1", value = 0)),
               column(2, numericInput("mgdgion2", "ion2", value = 0))),
             fluidRow(
               column(3, htmlOutput("mgdgsn1")),
               column(3, htmlOutput("mgdgsn2")),
               column(3, verbatimTextOutput("mgdgsum"))
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
             tags$li("[M + Na - sn2]+")),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - H]-"),
             br(),
             tags$li("[M - H - (sn1 - H2O)]-"),
             tags$li("[sn1 - H]-"),
             br(),
             tags$li("[M - H - (sn2 - H2O)]-"),
             tags$li("[sn2 - H]-")
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
             fluidRow(h4("ESI+"), 
                      column(2, strong("[M+NH4]+")),
                      column(6, verbatimTextOutput("dgdgfragpos"))),
             fluidRow(column(2, strong("[M+Na]+")), 
                      column(3, numericInput("dgdgion1pos", "ion1", value = 0)),
                      column(3, numericInput("dgdgion2pos", "ion2", value = 0))),
             fluidRow(
               column(2, ""),
               column(3, verbatimTextOutput("dgdgsn1pos")),
               column(3, verbatimTextOutput("dgdgsn2pos")),
               column(3, verbatimTextOutput("dgdgsumpos"))
             ),
             fluidRow(h4("ESI-"), verbatimTextOutput("dgdgfragneg")),
             fluidRow(
               column(2, numericInput("dgdgion1", "ion1", value = 0)),
               column(2, numericInput("dgdgion2", "ion2", value = 0))),
             fluidRow(
               column(3, htmlOutput("dgdgsn1")),
               column(3, htmlOutput("dgdgsn2")),
               column(3, verbatimTextOutput("dgdgsum"))
             )
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for DGDGs:"),
      column(3,
             strong("Positive [M+NH4]+:"),
             tags$li("[M + NH4 - 2(galactose)]+ + 0.984"),
             tags$li("[M + H - 2(galactose)]+"),
             br(),
             strong("Positive [M+Na]+:"),
             tags$li("[M + Na - sn1]+"),
             tags$li("[M + Na - galactose - sn1]+")
      ),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - H]-"),
             tags$li("[M - H - sn1]-")
      )
    ) # close DGDG
  ), #  close Glycerolypids
  
  ## Others ----
  navbarMenu(
    "Others",
    ### FFA ----
    tabPanel(
      "FFAs",
      h3("Free Fatty Acids (FFAs)"),
      column(4, strong("Formula"),
             fluidRow(
               column(2, numericInput("ffaC", "C", value = 18)),
               column(2, numericInput("ffadb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("ffaformula"))),
             fluidRow(),
             fluidRow(strong("m/z values"), verbatimTextOutput("ffamzvals1"))),
      column(1), 
      column(6, 
             h3("Commonly occuring product ions for Free Fatty Acids:"),
             p("STDs: there are only MS2 for [M-H]-, but any of them seems to be of good quality"),
             strong("Study samples:"),
             fluidRow(verbatimTextOutput("ffafrag"))
      ),
      fluidRow(),
      hr(),
      h3("FFAs reduced"),
      column(4, strong("Formula"),
             fluidRow(
               column(2, numericInput("ffarC", "C", value = 18)),
               column(2, numericInput("ffardb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("ffarformula"))),
             fluidRow(),
             fluidRow(strong("m/z values"), verbatimTextOutput("ffarmzvals1"))),
      column(1),
      column(6,
             verbatimTextOutput("ffarfrag")),
      fluidRow(),
      hr(),
      h3("FFAs + Phenylacetic acid"),
      column(4, strong("Formula"),
             fluidRow(
               column(2, numericInput("ffapaC", "C", value = 18)),
               column(2, numericInput("ffapadb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("ffapaformula"))),
             fluidRow(),
             fluidRow(strong("m/z values"), verbatimTextOutput("ffapamzvals1"))),
      column(1),
      column(6,
             verbatimTextOutput("ffapafrag")),
    ), # close tab FFA
    
    ### Carnitines ----
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
             tags$li("[M + H - trimethylamine]+"),
             tags$li("[sn + H - H2O]+")),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - H - carnitine]- / [sn - H]-")
      )
    ), # close tab Carnitines
    
    ### Ceramides ----
    tabPanel(
      "Ceramides",
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("cerC", "C", value = 15)),
               column(2, numericInput("cerdb", "db", value = 0))
             ),
             column(4, fluidRow(verbatimTextOutput("cerformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("cermzvals1")),
             fluidRow(verbatimTextOutput("cermzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("cerfragpos")),
             fluidRow(h4("ESI-"), verbatimTextOutput("cerfragneg"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for ceramides:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - H2O]+")),
      column(3, 
             strong("Negative [M-H+HCOOH]-:"),
             tags$li("[M - H]-")
      )
    ), # close tab Ceramides
    
    ### Sphingomyelines ----
    tabPanel(
      "Sphingomyelines",
      column(4, h3("Formula"),
             fluidRow(
               column(2, numericInput("smC", "C", value = 36)),
               column(2, numericInput("smdb", "db", value = 2))
             ),
             column(4, fluidRow(verbatimTextOutput("smformula"))),
             fluidRow(),
             fluidRow(h3("m/z values"), verbatimTextOutput("smmzvals1")),
             fluidRow(verbatimTextOutput("smmzvals2"))
      ),
      column(1), 
      column(6, h3("MS2"),
             fluidRow(h4("ESI+"), verbatimTextOutput("smfragpos"))
      ),
      fluidRow(),
      hr(),
      h3("Commonly occuring product ions for Sphingomyelines:"),
      column(3,
             strong("Positive [M+H]+:"),
             tags$li("[M + H - H2O]+"))
    ) # close tab Sphingomyelines
  ), # close Others
  
  
  ## MS2 library ----
  navbarMenu("MS2 in-house library",
             tabPanel("MS/MS spectras",
                      sidebarLayout(
                        sidebarPanel(
                          sliderInput(
                            "mz", "m/z zoom:", min = 50, max = 2000, 
                            value = c(50, 1000), step = 10),
                          sliderInput(
                            "int", "intensity threshold:", min = 0, max = 100, 
                            value = 10, step = 1),
                          fluidRow(
                            strong("Filter MS2 spectra that contain the fragment:")
                          ),
                          fluidRow(
                            column(5, numericInput(
                              "frag", "With an m/z value", 
                              value = 0)),
                            column(2, numericInput("ppm", "ppm", value = 10)),
                            column(5, numericInput(
                              "intx", "With an intensity >", value = 1))
                          ),
                          fluidRow(
                            radioButtons(
                              "radio", label = "", 
                              choices = list("Fragment" = 1, "Neutral loss" = 2),
                              selected = 1)
                          )
                        ),
                        mainPanel(
                          fluidRow(dataTableOutput("table")),
                          fluidRow(
                            column(6, plotOutput("plot_MS2")),
                            column(6, plotOutput("plot_NL"))
                          ))
                      )),
             tabPanel("Correlations",
                      sidebarLayout(
                        sidebarPanel(
                          fluidRow(fileInput("file1", "Choose TXT file")),
                          fluidRow(numericInput("precursor", "Precursor m/z", 
                                                value = 517.42540)),
                          fluidRow(checkboxGroupInput(
                            "pol", "Polarity",
                            choices = list("POS" = 1, "NEG" = 0),
                            selected = 0)),
                          fluidRow(plotOutput("ms2_x")),
                          fluidRow(actionButton("go", "Go"))),
                        mainPanel(
                          fluidRow(
                            column(6, dataTableOutput("spectra")),
                            column(6, plotOutput("ms2_spectra"))
                          ),
                          fluidRow(
                            column(6, dataTableOutput("spectra_nl")),
                            column(6, plotOutput("nl_spectra"))
                          )
                        )
                      )),
             tabPanel(
               "MS/MS spectras (STDs)",
               sidebarLayout(
                 sidebarPanel(
                   #sliderInput(
                   #   "mz", "m/z zoom:", min = 50, max = 2000, 
                   #   value = c(50, 1000), step = 10),
                   sliderInput(
                     "int_std", "intensity threshold:", min = 0, max = 100, 
                     value = 10, step = 1),
                   fluidRow(
                     strong("Filter MS2 spectra that contain the fragment:")
                   ),
                   fluidRow(
                     column(5, numericInput(
                       "frag_std", "With an m/z value", 
                       value = 0)),
                     column(2, numericInput("ppm_std", "ppm", value = 10)),
                     column(5, numericInput(
                       "intx_std", "With an intensity >", value = 1))
                   ),
                   fluidRow(
                     radioButtons(
                       "radio_std", label = "", 
                       choices = list("Fragment" = 1, "Neutral loss" = 2),
                       selected = 1)
                   )
                 ),
                 mainPanel(
                   fluidRow(dataTableOutput("table_std")),
                   fluidRow(
                     column(6, plotOutput("plot_MS2_std")),
                     column(6, plotOutput("plot_NL_std"))
                   ))
               )),
  ),
  
  ## RT prediction ----
  tabPanel("Compound prediction from RT",
           numericInput("rt", "RT", value = 15),
           numericInput("mzs", "m/z", value = 545.455),
           hr(),
           column(6, verbatimTextOutput("cmps_rts")),
           column(6, verbatimTextOutput("cmps_mzs"))
           
  ),
  
  ## Theoretical MS2 ----
  tabPanel("Theoretical MS2",
           sidebarLayout(
             sidebarPanel(
               selectInput("class", "Lipid class:",
                           choices = list("FFA" = "FFA",
                                          "FFA_red" = "FFA_red",
                                          "FFA_PA" = "FFA_PA",
                                          "CAR" = "CAR",
                                          "CER" = "CER",
                                          "CER_FA" = "CER_FA",
                                          "CER_hex" = "CER_hex",
                                          "CER_hex_FA" = "CER_hex_FA",
                                          "CER_lact" = "CER_lact",
                                          "SM" = "SM",
                                          "PA" = "PA",
                                          "mPA" = "mPA",
                                          "dmPA" = "dmPA",
                                          "LPC" = "LPC",
                                          "LPC_MS3" = "LPC_MS3",
                                          "PC" = "PC",
                                          "LPE" = "LPE",
                                          "PE" = "PE",
                                          "PG" = "PG",
                                          "LPI" = "LPI",
                                          "PI" = "PI",
                                          "PS" = "PS",
                                          "MGDG" = "MGDG",
                                          "MGDG_Na" = "MGDG_Na",
                                          "DGDG" = "DGDG",
                                          "DGDG_Na" = "DGDG_Na",
                                          "DAG" = "DAG",
                                          "TAG" = "TAG"),
                           selected = "PA"),
               fluidRow(
                 column(3, numericInput("C", "C", value = 36)),
                 column(3, numericInput("db", "db", value = 2))
               ),
               fluidRow(
                 column(4, selectInput("sn1", "sn1",
                                       choices = sn$sn,
                                       selected = "18:0")),
                 column(4, selectInput("sn2", "sn2",
                                       choices = sn$sn,
                                       selected = "18:2")),
                 column(4, selectInput("sn3", "sn3",
                                       choices = sn$sn,
                                       selected = "18:0"))
               )
             ),
             mainPanel(
               fluidRow(
                 column(6, plotOutput("thr_MS2_N")),
                 column(6, plotOutput("thr_MS2_P"))),
               fluidRow(
                 column(6, verbatimTextOutput("thr_MS2_tb_N")),
                 column(6, verbatimTextOutput("thr_MS2_tb_P")))
             )
           )
  )
  
  
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
    idx1 <- unlist(matchWithPpm(
      mass2mz(pamass(), "[M-H]-") - input$paion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pamass(), "[M-H]-") - input$paion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$paion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2, idx3)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pamass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pasn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pamass(), "[M-H]-") - input$paion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pamass(), "[M-H]-") - input$paion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$paion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2, idx3)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pamass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pasum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pamass(), "[M-H]-") - input$paion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pamass(), "[M-H]-") - input$paion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$paion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxa <- c(idx1, idx2, idx3)
    idx1 <- unlist(matchWithPpm(
      mass2mz(pamass(), "[M-H]-") - input$paion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pamass(), "[M-H]-") - input$paion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$paion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxb <- c(idx1, idx2, idx3)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
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
    idx1 <- unlist(matchWithPpm(
      mass2mz(mpamass(), "[M-H]-") - input$mpaion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mpaion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(mpamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$mpasn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(mpamass(), "[M-H]-") - input$mpaion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mpaion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(mpamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$mpasum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(mpamass(), "[M-H]-") - input$mpaion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mpaion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxa <- c(idx1, idx2)
    idx1 <- unlist(matchWithPpm(
      mass2mz(mpamass(), "[M-H]-") - input$mpaion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mpaion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxb <- c(idx1, idx2)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
  })
  
  ## PA dimethylated -----
  
  dmpafml <- reactive({
    paste0("C", input$dmpaC + 5, "H", 
           input$dmpaC*2 - (2 + 2*input$dmpadb) + 11, "O8P")
  })
  
  dmpamass <- reactive({
    MonoisotopicMass(formula = ListFormula(dmpafml()))
  })
  
  output$dmpaformula <- renderPrint({dmpafml()})
  
  output$dmpamzvals1 <- renderPrint({
    mass2mz(dmpamass(), adduct = c("[M+NH4]+", "[M-H]-"))
  })
  
  output$dmpamzvals2 <- renderPrint({
    mass2mz(dmpamass(), adduct = c("[2M+H]+", "[2M+NH4]+", "[2M-H]-"))
  })
  
  output$dmpafragpos <- renderPrint({
    tmp <- mass2mz(dmpamass(), "[M+H]+")
    tmp2 <- colnames(tmp)
    tmp <- c(as.numeric(unlist(mass2mz(dmpamass(), "[M+H]+"))) - 
               MonoisotopicMass(formula = ListFormula("H3PO4C2H4")),
             as.numeric(unlist(mass2mz(dmpamass(), "[M+H]+"))) - 
               MonoisotopicMass(formula = ListFormula("H3PO4")), tmp)
    names(tmp) <- c("[M+H-dmPA]+", "[M+H-PA]+", tmp2)
    tmp
  })
  
  output$dmpasn1 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(dmpamass(), "[M-H]-") - input$dmpaion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(dmpamass(), "[M-H]-") - input$dmpaion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$dmpaion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2, idx3)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(dmpamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(dmpamass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$dmpasn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(dmpamass(), "[M-H]-") - input$dmpaion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(dmpamass(), "[M-H]-") - input$dmpaion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$dmpaion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2, idx3)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(dmpamass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(dmpamass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$dmpasum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(dmpamass(), "[M-H]-") - input$dmpaion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(dmpamass(), "[M-H]-") - input$dmpaion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$dmpaion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxa <- c(idx1, idx2, idx3)
    idx1 <- unlist(matchWithPpm(
      mass2mz(dmpamass(), "[M-H]-") - input$dmpaion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(dmpamass(), "[M-H]-") - input$dmpaion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$dmpaion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxb <- c(idx1, idx2, idx3)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
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
  
  output$pgfragpos <- renderPrint({
    sprintf("%.5f", as.numeric(mass2mz(pgmass(),  adduct = c("[M+NH4]+"))) - 141.04)
  })
  
  output$pgsn1 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$pgion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion1 - MonoisotopicMass(formula = ListFormula("C3H8O3")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx <- c(idx1, idx2, idx3, idx4)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) - MonoisotopicMass(formula = ListFormula("C3H8O3"))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pgsn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$pgion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion2 - MonoisotopicMass(formula = ListFormula("C3H8O3")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx <- c(idx1, idx2, idx3, idx4)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) - MonoisotopicMass(formula = ListFormula("C3H8O3"))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pgsum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$pgion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion1 - MonoisotopicMass(formula = ListFormula("C3H8O3")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idxa <- c(idx1, idx2, idx3, idx4)
    idx1 <- unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pgmass(), "[M-H]-") - input$pgion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$pgion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pgmass(), "[M-H]-") - input$pgion2 - MonoisotopicMass(formula = ListFormula("C3H8O3")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idxb <- c(idx1, idx2, idx3, idx4)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
  })
  
  ## Lyso-PE ----
  lpefml <- reactive({
    paste0("C", input$lpeC + 5, "H", 
           input$lpeC*2 - (2 + 2*input$lpedb) + 14, "NO7P")
  })
  
  lpemass <- reactive({
    MonoisotopicMass(formula = ListFormula(lpefml()))
  })
  
  output$lpeformula <- renderPrint({lpefml()})
  
  output$lpemzvals <- renderPrint({
    mass2mz(lpemass(), adduct = c("[M+H]+", "[M-H]-"))
  })
  
  output$lpefragpos <- renderPrint({
    tmp <- c(sprintf("%.5f", as.numeric(mass2mz(lpemass(), "[M+H-H2O]+"))), 
             sprintf("%.5f", as.numeric(mass2mz(lpemass(), "[M+H]+")) - 
                       MonoisotopicMass(formula = ListFormula("C2H8NO4P"))))
    names(tmp) <- c("[M+H-H2O]+", "[M+H-PE]+")
    tmp
  })
  
  output$lpesn <- renderPrint({
    idx <- unlist(matchWithPpm(input$lpeion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    HTML(paste(sn$sn[idx],
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")),
               sep = '<br/>'))
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
  
  output$pemzvals2 <- renderPrint({
    c( mass2mz(pemass(), "[2M+H]+"), mass2mz(pemass(), "[2M-H]-"))
  })
  
  output$pefragpos <- renderPrint({
    round(as.numeric(mass2mz(pemass(), "[M+H]+")) - 
            MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 5)
  })
  
  output$pesn1 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pemass(), "[M-H]-") - input$peion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$peion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pemass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pesn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pemass(), "[M-H]-") - input$peion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$peion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pemass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pesum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pemass(), "[M-H]-") - input$peion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$peion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxa <- c(idx1, idx2)
    idx1 <- unlist(matchWithPpm(
      mass2mz(pemass(), "[M-H]-") - input$peion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$peion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxb <- c(idx1, idx2)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
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
  
  output$lpcfragpos <- renderPrint({
    tmp <- c(round(MonoisotopicMass(formula = ListFormula("C5H15NO4P")), 5),
             round(mass2mz(lpcmass(), "[M+H-H2O]+"), 5))
    names(tmp) <- c("[Phosphocholine]+", "[M+H-H2O]+")
    tmp
  })
  
  output$lpcfragneg <- renderPrint({
    tmp <- c(round(mass2mz(lpcmass(), "[M-H]-") - MonoisotopicMass(formula = ListFormula("CH2")), 5),
             mass2mz(sn$mass[sn$sn == paste0(input$lpcC, ":", input$lpcdb)], "[M-H]-"))
    names(tmp) <- c("[M-CH3]-", "[sn-H]-")
    tmp
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
  
  output$pcfragposna <- renderPrint({
    c(as.numeric(mass2mz(pcmass(), "[M+Na]+")) - MonoisotopicMass(formula = ListFormula("C3H9N")),
      as.numeric(mass2mz(pcmass(), "[M+Na]+")) - MonoisotopicMass(formula = ListFormula("C5H14NO4P")))
  })
  
  output$pcfragposh <- renderPrint({
    if(input$pcC < 36){
      tmp <- as.numeric(mass2mz(pcmass(), "[M+H]+") - 157.05)
      names(tmp) <- "[M+H-157.05]+"
      tmp
    } else {
      tmp <- as.numeric(mass2mz(pcmass(), "[M+H]+") - 183.0684)
      names(tmp) <- "[M+H-183.0684]+"
      tmp
    }
  })
  
  output$pcsn1 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pcmass(), "[M+H]+") - input$pcion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pcmass(), "[M+H]+") - input$pcion1, sn$mass, ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pcmass(), "[M+H]+") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pcmass(), "[M+H]+") - sn$mass[idx]),
               sep = '<br/>'))
  })
  
  output$pcsn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pcmass(), "[M+H]+") - input$pcion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pcmass(), "[M+H]+") - input$pcion2, sn$mass, ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pcmass(), "[M+H]+") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pcmass(), "[M+H]+") - sn$mass[idx]),
               sep = '<br/>'))
  })
  
  output$pcsum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pcmass(), "[M+H]+") - input$pcion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pcmass(), "[M+H]+") - input$pcion1, sn$mass, ppm = 10))
    idxa <- c(idx1, idx2)
    idx1 <- unlist(matchWithPpm(
      mass2mz(pcmass(), "[M+H]+") - input$pcion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pcmass(), "[M+H]+") - input$pcion2, sn$mass, ppm = 10))
    idxb <- c(idx1, idx2)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
  })
  
  output$pcfragneg <- renderPrint({
    tmp <- as.numeric(mass2mz(pcmass(), "[M-H]-")) - MonoisotopicMass(
      formula = ListFormula("CH2"))
    names(tmp) <- "[M-CH3]-"
    tmp
  })
  
  ## Lyso-PS ----
  lpsfml <- reactive({
    paste0("C", input$lpsC + 6, "H", 
           input$lpsC*2 - (2 + 2*input$lpsdb) + 14, "NO9P")
  })
  
  lpsmass <- reactive({
    MonoisotopicMass(formula = ListFormula(lpsfml()))
  })
  
  output$lpsformula <- renderPrint({lpsfml()})
  
  output$lpsmzvals1 <- renderPrint({
    mass2mz(lpsmass(), adduct = c("[M+H]+", "[M-H]-"))
  })
  
  output$lpsfragpos <- renderPrint({
    round(mass2mz(lpsmass(), "[M+H-H2O]+"), 5)
  })
  
  output$lpsfragneg <- renderPrint({
    as.numeric(round(mass2mz(lpsmass(), "[M-H]-") - 
                       MonoisotopicMass(formula = ListFormula("C3H5NO2")), 5))
  })
  
  ## PS ----
  psfml <- reactive({
    paste0("C", input$psC + 6, "H", 
           input$psC*2 - (2 + 2*input$psdb) + 12, "NO10P")
  })
  
  psmass <- reactive({
    MonoisotopicMass(formula = ListFormula(psfml()))
  })
  
  output$psformula <- renderPrint({psfml()})
  
  output$psmzvals <- renderPrint({
    mass2mz(psmass(), adduct = c("[M+H]+", "[M-H]-"))
  })
  
  output$psfragpos <- renderPrint({
    as.numeric(unlist(mass2mz(psmass(), "[M+H]+"))) - 
      MonoisotopicMass(formula = ListFormula("C3H8NO6P"))
  })
  
  output$psfragneg <- renderPrint({
    as.numeric(unlist(mass2mz(psmass(), "[M-H]-"))) - 
      MonoisotopicMass(formula = ListFormula("C3H5NO2"))
  })
  
  output$pssn1 <- renderPrint({
    idx <- unlist(matchWithPpm(
      input$psion1, 
      as.numeric(mass2mz(psmass(), "[M-H]-") - 
                   MonoisotopicMass(formula = ListFormula("C3H5NO2"))) - 
        sn$mass, ppm = 10))
    HTML(sn$sn[idx])
  })
  
  output$pssn2 <- renderPrint({
    idx <- unlist(matchWithPpm(
      input$psion1, 
      as.numeric(mass2mz(psmass(), "[M-H]-") - 
                   MonoisotopicMass(formula = ListFormula("C3H5NO2"))) - 
        sn$mass, ppm = 10))
    paste0(input$psC - sn$C[idx], ":", input$psdb - sn$db[idx])
  })
  
  ## Lyso-PI ----
  lpifml <- reactive({
    paste0("C", input$lpiC + 9, 
           "H", input$lpiC*2 - (2*input$lpidb) + 17, 
           "O12P")
  })
  
  lpimass <- reactive({
    MonoisotopicMass(formula = ListFormula(lpifml()))
  })
  
  output$lpiformula <- renderPrint({lpifml()})
  
  output$lpimzvals1 <- renderPrint({
    mass2mz(lpimass(), adduct = c("[M+H]+", "[M+NH4]+", "[M-H]-"))
  })
  
  #output$lpifragpos <- renderPrint({
  #  round(mass2mz(lpimass(), "[M+H-H2O]+"), 5)
  #})
  
  #output$lpifragneg <- renderPrint({
  #  as.numeric(round(mass2mz(lpimass(), "[M-H]-") - 
  #                     MonoisotopicMass(formula = ListFormula("C3H5NO2")), 5))
  #})
  
  
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
      adduct = c("[2M+NH4]+", "[2M-H]-")))
    tmp2 <- colnames(tmp)
    tmp <- c(as.numeric(unlist(mass2mz(dagmass(), "[M+H]+"))) + 
               MonoisotopicMass(formula = ListFormula("C2H7N")), tmp)
    names(tmp) <- c("[M+C2H8N]+", tmp2)
    #tmp <- tmp[1, 3, 2, 4]
    tmp
  })
  
  output$pifragpos <- renderPrint({
    tmp <- unlist(mass2mz(
      pimass(), 
      adduct = c("[M+H]+")))
    tmp2 <- colnames(tmp)
    tmp <- c(as.numeric(unlist(mass2mz(pimass(), "[M+NH4]+"))) - 
               MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984, tmp)
    names(tmp) <- c("[M+NH4-PI]+", tmp2)
    tmp <- tmp[c(2, 1)]
    tmp
  })
  
  output$pisn1 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pimass(), "[M-H]-") - input$piion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$piion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion1 - MonoisotopicMass(formula = ListFormula("C6H12O6")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx <- c(idx1, idx2, idx3, idx4)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) - MonoisotopicMass(formula = ListFormula("C6H12O6"))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pisn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pimass(), "[M-H]-") - input$piion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$piion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion2 - MonoisotopicMass(formula = ListFormula("C6H12O6")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx <- c(idx1, idx2, idx3, idx4)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - sn$mass[idx]),
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(pimass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O"))) - MonoisotopicMass(formula = ListFormula("C6H12O6"))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$pisum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pimass(), "[M-H]-") - input$piion1, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$piion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion1 - MonoisotopicMass(formula = ListFormula("C6H12O6")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idxa <- c(idx1, idx2, idx3, idx4)
    idx1 <- unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(mass2mz(pimass(), "[M-H]-") - input$piion2, sn$mass, ppm = 10))
    idx3 <- unlist(matchWithPpm(input$piion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx4 <-  unlist(matchWithPpm(
      mass2mz(pimass(), "[M-H]-") - input$piion2 - MonoisotopicMass(formula = ListFormula("C6H12O6")), 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idxb <- c(idx1, idx2, idx3, idx4)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
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
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$tagsn1],
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$tagsn2],
      unlist(mass2mz(tagmass(), "[M+H]+")) - sn$mass[sn$sn == input$tagsn3])
    tmp <- unique(tmp)
    names(tmp) <- unique(paste0(
      "[M+H-C", c(input$tagsn1, input$tagsn2, input$tagsn3), "]+"))
    tmp
  })
  
  output$tagformula <- renderPrint({tagfml()})
  
  output$tagmzvals1 <- renderPrint({
    mass2mz(tagmass(), "[M+NH4]+")
  })
  
  output$tagmzvals2 <- renderPrint({
    mass2mz(tagmass(), adduct = c("[M+H]+", "[2M+NH4]+"))
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
    idx1 <- unlist(matchWithPpm(
      mass2mz(mgdgmass(), "[M-H]-") - input$mgdgion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mgdgion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(mgdgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$mgdgsn2 <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(mgdgmass(), "[M-H]-") - input$mgdgion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mgdgion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idx <- c(idx1, idx2)
    HTML(paste(sn$sn[idx], 
               sprintf("%.5f", mass2mz(mgdgmass(), "[M-H]-") - (sn$mass[idx] - MonoisotopicMass(formula = ListFormula("H2O")))),
               sprintf("%.5f", mass2mz(sn$mass[idx], "[M-H]-")), 
               sep = '<br/>'))
  })
  
  output$mgdgsum <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      mass2mz(mgdgmass(), "[M-H]-") - input$mgdgion1, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mgdgion1, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxa <- c(idx1, idx2)
    idx1 <- unlist(matchWithPpm(
      mass2mz(mgdgmass(), "[M-H]-") - input$mgdgion2, 
      c(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))), ppm = 10))
    idx2 <- unlist(matchWithPpm(input$mgdgion2, mass2mz(sn$mass, "[M-H]-"), ppm = 10))
    idxb <- c(idx1, idx2)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
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
  
  output$dgdgsn1pos <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - input$dgdgion1pos, 
      sn$mass, ppm = 10))
    idx2 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - 
        MonoisotopicMass(formula = ListFormula("C6H10O5")) - input$dgdgion1pos, 
      sn$mass, ppm = 10))
    idx <- c(idx1, idx2)
    HTML(sn$sn[idx])
  })
  
  output$dgdgsn2pos <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - input$dgdgion2pos, 
      sn$mass, ppm = 10))
    idx2 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - 
        MonoisotopicMass(formula = ListFormula("C6H10O5")) - input$dgdgion2pos, 
      sn$mass, ppm = 10))
    idx <- c(idx1, idx2)
    HTML(sn$sn[idx])
  })
  
  output$dgdgsumpos <- renderPrint({
    idx1 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - input$dgdgion1pos, 
      sn$mass, ppm = 10))
    idx2 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - 
        MonoisotopicMass(formula = ListFormula("C6H10O5")) - input$dgdgion1pos, 
      sn$mass, ppm = 10))
    idxa <- c(idx1, idx2)
    idx1 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - input$dgdgion2pos, 
      sn$mass, ppm = 10))
    idx2 <- unlist(matchWithPpm(
      unlist(mass2mz(dgdgmass(), adduct = c("[M+Na]+"))) - 
        MonoisotopicMass(formula = ListFormula("C6H10O5")) - input$dgdgion2pos, 
      sn$mass, ppm = 10))
    idxb <- c(idx1, idx2)
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
  })
  
  output$dgdgfragneg <- renderPrint({
    round(mass2mz(dgdgmass(), "[M-H]-"), 5)
  })
  
  output$dgdgsn1 <- renderPrint({
    idx <- unlist(matchWithPpm(
      mass2mz(dgdgmass(), "[M-H]-") - input$dgdgion1, sn$mass, ppm = 10))
    HTML(sn$sn[idx])
  })
  
  output$dgdgsn2 <- renderPrint({
    idx <- unlist(matchWithPpm(
      mass2mz(dgdgmass(), "[M-H]-") - input$dgdgion2, sn$mass, ppm = 10))
    HTML(sn$sn[idx])
  })
  
  output$dgdgsum <- renderPrint({
    idxa <- unlist(matchWithPpm(
      mass2mz(dgdgmass(), "[M-H]-") - input$dgdgion1, sn$mass, ppm = 10))
    idxb <- unlist(matchWithPpm(
      mass2mz(dgdgmass(), "[M-H]-") - input$dgdgion2, sn$mass, ppm = 10))
    paste0(sn$C[idxa] + sn$C[idxb], ":", sn$db[idxa] + sn$db[idxb])
  })
  
  ## FFA ----
  ffafml <- reactive({
    paste0("C", input$ffaC, "H", 
           input$ffaC*2 - (2*input$ffadb), "O2")
  })
  
  ffamass <- reactive({
    MonoisotopicMass(formula = ListFormula(ffafml()))
  })
  
  output$ffaformula <- renderPrint({ffafml()})
  
  output$ffamzvals1 <- renderPrint({
    mass2mz(ffamass(), adduct = c("[M+H]+", "[M-H]-", "[M+CHO2]-"))
  })
  
  output$ffafrag <- renderPrint({
    mz <- as.numeric(mass2mz(ffamass(), "[M+CHO2]-") - MonoisotopicMass(formula = ListFormula("CO2")))
    names(mz) <- "[M-H+HCOOH-CO2]-"
    mz
  })
  
  
  ffarfml <- reactive({
    paste0("C", input$ffarC, 
           "H", input$ffarC*2 - (2*input$ffardb) + 2, 
           "O")
  })
  
  ffarmass <- reactive({
    MonoisotopicMass(formula = ListFormula(ffarfml()))
  })
  
  output$ffarformula <- renderPrint({ffarfml()})
  
  output$ffarmzvals1 <- renderPrint({
    mass2mz(ffarmass(), adduct = c("[M+H]+", "[M-H]-", "[M+CHO2]-"))
  })
  
  output$ffarfrag <- renderPrint({
    mass2mz(ffarmass(), "[M-H]-")
  })
  
  
  
  ffapafml <- reactive({
    paste0("C", input$ffapaC + 8, 
           "H", input$ffapaC*2 - (2*input$ffardb) + 8, 
           "O", 2 + 2)
  })
  
  ffapamass <- reactive({
    MonoisotopicMass(formula = ListFormula(ffapafml()))
  })
  
  output$ffapaformula <- renderPrint({ffapafml()})
  
  output$ffapamzvals1 <- renderPrint({
    mass2mz(ffapamass(), adduct = c("[M+H]+", "[M-H]-", "[M+CHO2]-"))
  })
  
  output$ffapafrag <- renderPrint({
    mz <- as.numeric(mass2mz(ffapamass(), "[M-H]-") - MonoisotopicMass(formula = ListFormula("C8H8O2")))
    names(mz) <- "[M-H-PA]-"
    mz
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
    tmp <- c(round(as.numeric(mass2mz(carmass(), "[M+H]+")) - 
                     MonoisotopicMass(formula = ListFormula("C3H9N")), 5),
             round(as.numeric(mass2mz(
               MonoisotopicMass(formula = ListFormula(
                 paste0("C", input$carC, "H", 
                        input$carC*2 - (2*input$cardb), "O2"))), "[M+H-H2O]+")), 
               5))
    names(tmp) <- c("[M+H-TMA]+", "[sn+H-H2O]+")
    tmp
  })
  
  output$carfragneg <- renderPrint({
    round(as.numeric(mass2mz(carmass(), "[M-H]-")) - 
            MonoisotopicMass(formula = ListFormula("C7H13NO2")), 5)
  })
  
  ## Ceramides ----
  cerfml <- reactive({
    paste0("C", input$cerC, "H", 
           input$cerC*2 - (2*input$cerdb) + 1, "NO3")
  })
  
  cermass <- reactive({
    MonoisotopicMass(formula = ListFormula(cerfml()))
  })
  
  output$cerformula <- renderPrint({cerfml()})
  
  output$cermzvals1 <- renderPrint({
    mass2mz(cermass(), adduct = c("[M+H]+", "[M+CHO2]-"))
  })
  
  output$cerfragpos <- renderPrint({
    mass2mz(cermass(), "[M+H-H2O]+")
  })
  
  output$cerfragneg <- renderPrint({
    mass2mz(cermass(), "[M-H]-")
  })
  
  ## Sphingomyelines ----
  smfml <- reactive({
    paste0("C", input$smC + 4, 
           "H", input$smC*2 - (2*input$smdb) + 11, 
           "N2O6P")
  })
  
  smmass <- reactive({
    MonoisotopicMass(formula = ListFormula(smfml()))
  })
  
  output$smformula <- renderPrint({smfml()})
  
  output$smmzvals1 <- renderPrint({
    mass2mz(smmass(), adduct = c("[M+H]+", "[M+CHO2]-"))
  })
  
  output$smfragpos <- renderPrint({
    mass2mz(smmass(), "[M+H-H2O]+")
  })
  
  # MS2 library ----
  dbx <- reactive({
    if(input$frag > 0){
      if(input$radio == 1){
        tmp <- intensity(sps_ms2)[matchWithPpm(
          input$frag, mz(sps_ms2), ppm = input$ppm)[[1]]]/max(intensity(sps_ms2))
        tmp <- lapply(tmp, function(x) if (length(x) == 0) {0} else {x})
        sps_ms2 <- sps_ms2[unlist(tmp>(input$intx/100))]
      } else if(input$radio == 2){
        #tmp <- unlist(matchWithPpm(input$frag, abs(mz(sps_ms2_nl)), ppm = input$ppm))[[1]]
        tmp <- intensity(sps_ms2)[matchWithPpm(
          input$frag, abs(mz(sps_ms2_nl)), ppm = input$ppm)[[1]]]/max(intensity(sps_ms2))
        tmp <- lapply(tmp, function(x) if (length(x) == 0) {0} else {x})
        #sps_ms2 <- sps_ms2[which(tmp > 0)]
        sps_ms2 <- sps_ms2[unlist(tmp>(input$intx/100))]
      }
    }
    
    db <- data.frame(cbind(class = sps_ms2$class, 
                           compound = sps_ms2$name, 
                           adduct = sps_ms2$adduct))
    db$class[is.na(db$class)] <- "unknown"
    db$precursor <- as.numeric(sprintf("%.5f", precursorMz(sps_ms2)))
    db <- db[order(db$class, db$compound, db$adduct),]
    db
  })
  
  output$table <- renderDataTable(datatable({
    dbx()
  }, rownames = FALSE))
  
  output$plot_MS2 <- renderPlot({
    db <- dbx()
    i <- input$table_rows_selected
    if(length(i) == 1){
      j <- which(sps_ms2$name == db$compound[i] & 
                   sps_ms2$adduct == db$adduct[i])
      plot(unlist(mz(sps_ms2[j])),
           unlist(intensity(sps_ms2[j])), type = "h", 
           xlab = "m/z", ylab = "relative intensity", 
           xlim = c(input$mz[1],#min(unlist(mz(sps_ms2[j])))-10, 
                    input$mz[2]), #max(unlist(mz(sps_ms2[j])))+10), 
           ylim = c(0, 105),
           main = paste("MS2 for m/z", 
                        sprintf("%.5f", precursorMz(sps_ms2[j])), "@", 
                        sprintf("%.2f", rtime(sps_ms2[j])/60), "min"))
      idx <- unlist(intensity(sps_ms2[j])) > input$int
      text(unlist(mz(sps_ms2[j]))[idx],
           unlist(intensity(sps_ms2[j]))[idx], 
           sprintf("%.5f", unlist(mz(sps_ms2[j])))[idx], pos = 3 
           #offset = -1, pos = 2, srt = -30
      )
    }
  })
  
  output$plot_NL <- renderPlot({
    db <- dbx()
    i <- input$table_rows_selected
    if(length(i) == 1){
      j <- which(sps_ms2$name == db$compound[i] & 
                   sps_ms2$adduct == db$adduct[i])
      mzvals <- as.numeric(unlist(mz(sps_ms2[j]))) - precursorMz(sps_ms2[j])
      plot(mzvals, unlist(intensity(sps_ms2[j])), type = "h", 
           xlab = "m/z", ylab = "relative intensity", main = "Neutral Loss",
           xlim = c(50 - precursorMz(sps_ms2[j]), 
                    0),#min(mzvals) - 10, max(mzvals) + 10), 
           ylim = c(0, 105))
      idx <- unlist(intensity(sps_ms2[j])) > input$int
      text(mzvals[idx], unlist(intensity(sps_ms2[j]))[idx], 
           sprintf("%.5f", mzvals)[idx], pos = 3)
    }
  })
  
  x_spdx <- reactive({
    req(input$file1)
    df <- read.table(input$file1$datapath)
    df <- df[order(df[,1]),]
    x_spd <- DataFrame(
      msLevel = 2L,
      polarity = 1L,
      precursorMz = input$precursor,
      id = "x",
      name = "x")
    x_spd$mz <- list(df[,1])
    x_spd$intensity <- list(df[,2])
    x_spd <- Spectra(x_spd)
  })
  
  tbx <- eventReactive(input$go, {
    x_spd <- x_spdx()
    spdx <- filterPolarity(sps_ms2, input$pol)
    #spdx <- spdx[grep("Unk-00", spdx$name)]
    c_spd <- c(x_spd, spdx)
    tb <- cbind(c_spd$name,  Spectra::compareSpectra(c_spd, ppm = 50)[,1],
                c_spd$adduct)
    tb <- tb[-1,]
    colnames(tb) <- c("name", "corr", "adduct")
    tb[,2] <- sprintf("%.3f", round(as.numeric(tb[,2]), 3))
    return(tb)
  })
  
  #tbx <- reactive({
  #  x_spd <- x_spdx()
  #  spdx <- filterPolarity(sps_ms2, input$pol)
  #  spdx <- spdx[grep("Unk-00", spdx$name)]
  #  c_spd <- c(x_spd, spdx)
  #  tb <- cbind(c_spd$name,  Spectra::compareSpectra(c_spd, ppm = 50)[,1],
  #              c_spd$adduct)
  #  tb <- tb[-1,]
  #  colnames(tb) <- c("name", "corr", "adduct")
  #  tb[,2] <- sprintf("%.3f", round(as.numeric(tb[,2]), 3))
  #  return(tb)
  #})
  
  #tby <- reactive({
  tby <- eventReactive(input$go, {
    x_spd <- x_spdx()
    spdx <- filterPolarity(sps_ms2, input$pol)
    #spdx <- spdx[grep("Unk-00", spdx$name)]
    c_spd <- c(x_spd, spdx)
    c_spd <- applyProcessing(c_spd)
    mz(c_spd@backend) <- mz(c_spd) - precursorMz(c_spd)
    tb <- cbind(c_spd$name,  Spectra::compareSpectra(c_spd, tolerance = 0.01)[,1],
                c_spd$adduct)
    tb <- tb[-1,]
    colnames(tb) <- c("name", "corr", "adduct")
    tb[,2] <- sprintf("%.3f", round(as.numeric(tb[,2]), 3))
    return(tb) 
  })
  
  output$spectra <- renderDataTable(datatable({
    tb1 <- tbx()
    #tb <- tb[order(tb[,"corr"], decreasing = TRUE), ]
    return(tb1)
  }))
  
  output$spectra_nl <- renderDataTable(datatable({
    tb2 <- tby()
    #tb <- tb[order(tb[,"corr"], decreasing = TRUE), ]
    return(tb2)
  }))
  
  output$ms2_spectra <- renderPlot({
    x_spd <- x_spdx()
    tb <- tbx()
    
    i <- input$spectra_rows_selected
    if(length(i) == 1){
      j <- which(sps_ms2$name == tb[i, "name"] & 
                   sps_ms2$adduct == tb[i, "adduct"] &
                   sps_ms2$polarity == input$pol)
      plotSpectraMirror(x_spd, sps_ms2[j], tolerance = 0.01,
                        labels = label_fun, labelPos = 2, labelOffset = 0.2,
                        labelSrt = -30)
      grid()
    }
  })
  
  output$nl_spectra <- renderPlot({
    x_spd <- x_spdx()
    tb <- tby()
    x_spd <- applyProcessing(x_spd)
    mz(x_spd@backend) <- mz(x_spd) - precursorMz(x_spd)
    sps_ms2 <- applyProcessing(sps_ms2)
    mz(sps_ms2@backend) <- mz(sps_ms2) - precursorMz(sps_ms2)
    
    i <- input$spectra_nl_rows_selected
    if(length(i) == 1){
      j <- which(sps_ms2$name == tb[i, "name"] & 
                   sps_ms2$adduct == tb[i, "adduct"] &
                   sps_ms2$polarity == input$pol)
      plotSpectraMirror(x_spd, sps_ms2[j], tolerance = 0.01,
                        labels = label_fun, labelPos = 2, labelOffset = 0.2,
                        labelSrt = -30)
      grid()
    }
  })
  
  # MS2 library (STDs) ----
  dbx_std <- reactive({
    if(input$frag_std > 0){
      if(input$radio_std == 1){
        tmp <- intensity(sps_ms2_std)[matchWithPpm(
          input$frag_std, mz(sps_ms2_std), ppm = input$ppm_std)[[1]]]/max(intensity(sps_ms2_std))
        tmp <- lapply(tmp, function(x) if (length(x) == 0) {0} else {x})
        sps_ms2_std <- sps_ms2_std[unlist(tmp>(input$intx_std/100))]
      } else if(input$radio_std == 2){
        tmp <- intensity(sps_ms2_std)[matchWithPpm(
          input$frag_std, abs(mz(sps_ms2_std_nl)), ppm = input$ppm_std)[[1]]]/max(intensity(sps_ms2_std))
        tmp <- lapply(tmp, function(x) if (length(x) == 0) {0} else {x})
        sps_ms2_std <- sps_ms2_std[unlist(tmp>(input$intx_std/100))]
      }
    }
    
    db <- data.frame(cbind(class = sps_ms2_std$class, 
                           compound = sps_ms2_std$name, 
                           adduct = sps_ms2_std$adduct))
    db$class[is.na(db$class)] <- "unknown"
    db$precursor <- as.numeric(sprintf("%.5f", precursorMz(sps_ms2_std)))
    db <- db[order(db$class, db$compound, db$adduct),]
    db
  })
  
  output$table_std <- renderDataTable(datatable({
    dbx_std()
  }, rownames = FALSE))
  
  output$plot_MS2_std <- renderPlot({
    db <- dbx_std()
    k <- input$table_std_rows_selected
    if(length(k) == 1){
      j <- which(sps_ms2_std$name == db$compound[k] & 
                   sps_ms2_std$adduct == db$adduct[k])
      plot(unlist(mz(sps_ms2_std[j])),
           unlist(intensity(sps_ms2_std[j])), type = "h", 
           xlab = "m/z", ylab = "relative intensity", 
           xlim = c(50,#input$mz[1],#min(unlist(mz(sps_ms2_std[j])))-10, 
                    precursorMz(sps_ms2_std[j])),#input$mz[2]), #max(unlist(mz(sps_ms2_std[j])))+10), 
           ylim = c(0, 105),
           main = paste("MS2 for m/z", 
                        sprintf("%.5f", precursorMz(sps_ms2_std[j])), "@", 
                        sprintf("%.2f", rtime(sps_ms2_std[j])/60), "min"))
      idx <- unlist(intensity(sps_ms2_std[j])) > input$int_std
      text(unlist(mz(sps_ms2_std[j]))[idx],
           unlist(intensity(sps_ms2_std[j]))[idx], 
           sprintf("%.5f", unlist(mz(sps_ms2_std[j])))[idx], pos = 3 
           #offset = -1, pos = 2, srt = -30
      )
    }
  })
  
  output$plot_NL_std <- renderPlot({
    db <- dbx_std()
    k <- input$table_std_rows_selected
    if(length(k) == 1){
      j <- which(sps_ms2_std$name == db$compound[k] & 
                   sps_ms2_std$adduct == db$adduct[k])
      mzvals <- as.numeric(unlist(mz(sps_ms2_std[j]))) - precursorMz(sps_ms2_std[j])
      plot(mzvals, unlist(intensity(sps_ms2_std[j])), type = "h", 
           xlab = "m/z", ylab = "relative intensity", main = "Neutral Loss",
           xlim = c(50 - precursorMz(sps_ms2_std[j]), 
                    0),#min(mzvals) - 10, max(mzvals) + 10), 
           ylim = c(0, 105))
      idx <- unlist(intensity(sps_ms2_std[j])) > input$int_std
      text(mzvals[idx], unlist(intensity(sps_ms2_std[j]))[idx], 
           sprintf("%.5f", mzvals)[idx], pos = 3)
    }
  })
  
  ## RT ----
  output$cmps_rts <- renderPrint({
    rts <- as.data.frame(rbind(
      c("FFA", 8.469, 0.311, 0),
      c("PA", 2.592, 0.403, -0.643),
      c("mPA", 3.591, 0.365, -0.616),
      c("dmPA", 4.895, 0.325, -0.580),
      #c("LPC", 9.720, 0, -1.380),
      c("PC", 1.640, 0.382, -0.624),
      c("PE", 4.269, 0.349, -0.605),
      c("PG", 2.418, 0.364, -0.602),
      c("PI", -0.971, 0.414, -0.644),
      c("MGDG", 3.577, 0.341, -0.586),
      c("DGDG", -1.859, 0.402, -0.653),
      c("DAG", 6.878, 0.329, -0.539),
      c("TAG", 14.006, 0.151, -0.280)))
    colnames(rts) <- c("class", "intrs", "Cx", "dbx")
    rts$intrs <- as.numeric(rts$intrs)
    rts$Cx <- as.numeric(rts$Cx)
    rts$dbx <- as.numeric(rts$dbx)
    
    rt <- input$rt
    tb <- data.frame(matrix(nrow = 0, ncol = nrow(rts)))
    for(i in 0:10){
      tb <- rbind(tb, round((rt - rts$intrs - i*rts$dbx)/rts$Cx, 1))
    }
    colnames(tb) <- rts$class
    rownames(tb) <- 0:10
    tb[round((tb %% 1)*10) < 3] <- round(tb[round((tb %% 1)*10) < 3] )
    tb[round((tb %% 1)*10) > 7] <- round(tb[round((tb %% 1)*10) > 7] )
    tb[(tb %% 1)*10 > 0] <- NA
    tb <- tb[rowSums(is.na(tb)) != ncol(tb), ]
    for(i in seq(ncol(tb))){
      idx <- which(!is.na(tb[,i]))
      for(j in idx){
        if(colnames(tb)[i] == "FFA"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i])*2 - (2*as.numeric(rownames(tb)[j])), "O2")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("FFA", (tb[j,i]), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+H]+", "[M+CHO2]-"))))
          
        } else if(colnames(tb)[i] == "PA"){
          i.fml <- paste0("C", (tb[j,i]), 
                          "H", (tb[j,i] - 3)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 7, "O8P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("PA", (tb[j,i] - 3), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "mPA"){
          i.fml <- paste0("C", (tb[j,i]), 
                          "H", (tb[j,i] - 4)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 9, "O8P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("mPA", (tb[j,i] - 4), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "dmPA"){
          i.fml <- paste0("C", (tb[j,i]), 
                          "H", (tb[j,i] - 5)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 11, "O8P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("mPA", (tb[j,i] - 5), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "PC"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 8)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 18, "NO8P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("PC", (tb[j,i] - 8), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+H]+", "[M+CHO2]-"))))
          
        } else if(colnames(tb)[i] == "PE"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 5)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 12, "NO8P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("PE", (tb[j,i] - 5), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+H]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "PG"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 6)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 13, "O10P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("PG", (tb[j,i] - 6), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "PI"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 9)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 17, "O13P")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("PI", (tb[j,i] - 9), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M-H]-"))))
          
        } else if(colnames(tb)[i] == "MGDG"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 9)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 16, "O10")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("MGDG", (tb[j,i] - 9), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M+CHO2]-"))))
          
        } else if(colnames(tb)[i] == "DGDG"){
          i.fml <- paste0("C", (tb[j,i]), "H", (tb[j,i] - 15)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 26, "O15")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("DGDG", (tb[j,i] - 15), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M+CHO2]-"))))
          
        } else if(colnames(tb)[i] == "DAG"){
          i.fml <- paste0("C", (tb[j,i] - 3), "H", (tb[j,i] - 3)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 6, "O5")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("DAG", (tb[j,i] - 3), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M+CHO2]-"))))
          
        } else if(colnames(tb)[i] == "TAG"){
          i.fml <- paste0("C", (tb[j,i] - 3), "H", (tb[j,i] - 3)*2 - (2 + 2*as.numeric(rownames(tb)[j])) + 5, "O6")
          i.mass <- MonoisotopicMass(formula = ListFormula(i.fml))
          print(cbind(paste0("TAG", (tb[j,i] - 3), ":", rownames(tb)[j]), mass2mz(i.mass, adduct = c("[M+NH4]+", "[M+CHO2]-"))))
          
        }
      }
    }
    
  })
  
  output$cmps_mzs <- renderPrint({
    dt <- matrix(nrow = 0, ncol = 3)
    for(c in 10:30){
      for(db in 0:4){
        dt <- rbind(
          dt,
          c(paste0("CAR", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 7, 
                     "H", c*2 - (2*db) + 13, "NO4"))), "[M+CHO2]-")),
          c(paste0("CAR", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 7, 
                     "H", c*2 - (2*db) + 13, "NO4"))),"[M+H]+")),
          c(paste0("LPC", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 8, 
                     "H", c*2 - (2*db) + 18, "NO7P"))), "[M+CHO2]-")),
          c(paste0("LPC", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 8, 
                     "H", c*2 - (2*db) + 18, "NO7P"))),"[M+H]+")),
          c(paste0("LPE", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 12, "NO7P"))), "[M-H]-")),
          c(paste0("LPE", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 12, "NO7P"))),"[M+H]+")),
          c(paste0("LPG", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 13, "O9P"))), "[M-H]-")),
          c(paste0("LPG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 13, "O9P"))),"[M+NH4]+")),
          c(paste0("LPI", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 17, "O12P"))), "[M-H]-")),
          c(paste0("LPI", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 17, "O12P"))),"[M+NH4]+")),
          c(paste0("LPS", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 12, "NO9P"))), "[M-H]-")),
          c(paste0("LPS", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 12, "NO9P"))),"[M+H]+"))
        )
      }}
    for(c in 10:40){
      for(db in 0:4){
        dt <- rbind(
          dt, 
          c(paste0("FFA", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, "H", c*2 - (2*db), "O2"))), "[M+CHO2]-"))
        )
      }}
    for(c in 20:50){
      for(db in 0:6){
        dt <- rbind(
          dt, 
          c(paste0("PA", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 5, "O8P"))), "[M-H]-")),
          c(paste0("PA", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 5, "O8P"))),"[M+NH4]+")),
          c(paste0("mPA", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 4, 
                     "H", c*2 - (2*db) + 7, "O8P"))), "[M-H]-")),
          c(paste0("mPA", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 4, 
                     "H", c*2 - (2*db) + 7, "O8P"))),"[M+NH4]+")),
          c(paste0("dmPA", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 9, "O8P"))), "[M-H]-")),
          c(paste0("dmPA", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 9, "O8P"))),"[M+NH4]+")),
          c(paste0("PC", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 8, 
                     "H", c*2 - (2*db) + 16, "NO8P"))), "[M+CHO2]-")),
          c(paste0("PC", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 8, 
                     "H", c*2 - (2*db) + 16, "NO8P"))),"[M+H]+")),
          c(paste0("PE", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 10, "NO8P"))), "[M-H]-")),
          c(paste0("PE", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 10, "NO8P"))),"[M+H]+")),
          c(paste0("PG", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 11, "O10P"))), "[M-H]-")),
          c(paste0("PG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 11, "O10P"))),"[M+NH4]+")),
          c(paste0("PI", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 15, "O13P"))), "[M-H]-")),
          c(paste0("PI", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 15, "O13P"))),"[M+NH4]+")),
          c(paste0("PS", c, ":", db), 
            "[M-H]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 10, "NO10P"))), "[M-H]-")),
          c(paste0("PS", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 10, "NO10P"))),"[M+H]+")),
          c(paste0("SM", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 13, "N2O6P"))), "[M+CHO2]-")),
          c(paste0("SM", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 5, 
                     "H", c*2 - (2*db) + 13, "N2O6P"))),"[M+H]+")),
          c(paste0("CER", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO3"))), "[M+CHO2]-")),
          c(paste0("CER", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO3"))),"[M+H]+")),
          c(paste0("CER_gluc", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 11, "NO8"))), "[M+CHO2]-")),
          c(paste0("CER_gluc", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 11, "NO8"))),"[M+H]+")),
          c(paste0("CER_lact", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 12, 
                     "H", c*2 - (2*db) + 21, "NO13"))), "[M+CHO2]-")),
          c(paste0("CER_lact", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 12, 
                     "H", c*2 - (2*db) + 21, "NO13"))),"[M+H]+")),
          c(paste0("CER_hex3", c , ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 18, 
                     "H", c*2 - (2*db) + 31, "NO18"))), "[M+CHO2]-")),
          c(paste0("CER_hex3", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 18, 
                     "H", c*2 - (2*db) + 31, "NO18"))),"[M+H]+")),
          c(paste0("dihydroCER", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 3, "NO3"))), "[M+CHO2]-")),
          c(paste0("dihydroCER", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 3, "NO3"))),"[M+H]+")),
          c(paste0("dihydroCER_gluc", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 13, "NO8"))), "[M+CHO2]-")),
          c(paste0("dihydroCER_gluc", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 6, 
                     "H", c*2 - (2*db) + 13, "NO8"))),"[M+H]+")),
          c(paste0("dihydroCER_lact", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 12, 
                     "H", c*2 - (2*db) + 23, "NO13"))), "[M+CHO2]-")),
          c(paste0("dihydroCER_lact", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 12, 
                     "H", c*2 - (2*db) + 23, "NO13"))),"[M+H]+")),
          c(paste0("deoxyCER", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO2"))), "[M+CHO2]-")),
          c(paste0("deoxyCER", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO2"))),"[M+H]+")),
          c(paste0("phytoCER", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO4"))), "[M+CHO2]-")),
          c(paste0("phytoCER", c, ":", db), 
            "[M+H]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c, 
                     "H", c*2 - (2*db) + 1, "NO4"))),"[M+H]+")),
          c(paste0("MGDG", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 14, "O10"))), "[M+CHO2]-")),
          c(paste0("MGDG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 9, 
                     "H", c*2 - (2*db) + 14, "O10"))),"[M+NH4]+")),
          c(paste0("DGDG", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 15, 
                     "H", c*2 - (2*db) + 24, "O15"))), "[M+CHO2]-")),
          c(paste0("DGDG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 15, 
                     "H", c*2 - (2*db) + 24, "O15"))),"[M+NH4]+")),
          c(paste0("DAG", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 4, "O5"))), "[M+CHO2]-")),
          c(paste0("DAG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 4, "O5"))),"[M+NH4]+"))
        )
      }}
    for(c in 35:70){
      for(db in 0:10){
        dt <- rbind(
          dt, 
          c(paste0("TAG", c, ":", db), 
            "[M+CHO2]-",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 2, "O6"))), "[M+CHO2]-")),
          c(paste0("TAG", c, ":", db), 
            "[M+NH4]+",
            mass2mz(MonoisotopicMass(formula = ListFormula(
              paste0("C", c + 3, 
                     "H", c*2 - (2*db) + 2, "O6"))),"[M+NH4]+"))
        )
      }
    }
    dt <- data.frame(dt)
    colnames(dt) <- c("compound", "adduct", "mz")
    dt$mz <- as.numeric(dt$mz)
    xdt <- dt[unlist(matchWithPpm(input$mzs, dt$mz, ppm = 10)),]
    print(xdt)
  })
  
  ## Theoretical MS2 ----
  #### thr_MS2_NEG ----
  thr_MS2_NEG <- reactive({
    if(input$class == "FFA"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C - 1, "H", input$C*2 - (2*input$db), "HCOOH"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M-H-HCOOH-CO2]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "FFA_red"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db) + 2,
        "O"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M-H]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "FFA_PA"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db),
        "O2"))), "[M-H]-")
      dt <- data.frame(cbind(
        "mz" = as.numeric(i.mz), 
        "percent" = 100,
        "fragment" = paste0("[", input$C, ":", input$db, "-H]-")))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CAR"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db), "O2"))), "[M-H]-")
      dt <- data.frame(cbind(
        "mz" = as.numeric(i.mz), 
        "percent" = 100,
        "fragment" = paste0("[", input$C, ":", input$db, "-H]-")))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CER"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C - 1, 
        "H", input$C*2 - (2*input$db) + 1 - 4, 
        "N",
        "O", 3 - 1))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M-H-CH3OH]-"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db) + 1 - 2, 
        "N",
        "O", 3 - 1))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 15,
                              "fragment" = "[M-H-H2O]-"))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C - 1, 
        "H", input$C*2 - (2*input$db) + 1 - 4, 
        "N",
        "O", 3 - 2))), "[M-H]-")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 60,
                              "fragment" = "[M-H-CH4O2]-"))
      dt <- rbind(dt, dt3)
      
      i.mz4 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)) + 2, 
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 - (
          2*as.numeric(gsub(".*:", "", input$sn1))) + 3, 
        "N",
        "O", 2 - 1))), "[M-H]-")
      dt4 <- data.frame(cbind(
        "mz" = as.numeric(i.mz4), 
        "percent" = 80,
        "fragment" = paste0("[", input$sn1, "+C2H3N-O]-")))
      dt <- rbind(dt, dt4)
      
      i.mz5 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)) + 2, 
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 - (
          2*as.numeric(gsub(".*:", "", input$sn1))) + 3, 
        "N",
        "O", 2))), "[M-H]-")
      dt5 <- data.frame(cbind(
        "mz" = as.numeric(i.mz5), 
        "percent" = 30,
        "fragment" = paste0("[", input$sn1, "+C2H3N]-")))
      dt <- rbind(dt, dt5)
      
      i.mz6 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)), 
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 - 
          (2*as.numeric(gsub(".*:", "", input$sn1))), 
        "O", 2))), "[M-H]-")
      dt6 <- data.frame(cbind(
        "mz" = as.numeric(i.mz6), 
        "percent" = 40,
        "fragment" = paste0("[", input$sn1, "-H]-")))
      dt <- rbind(dt, dt6)
      
      i.mz6 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)), 
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 - 
          (2*as.numeric(gsub(".*:", "", input$sn1))) - 2, 
        "O", 2 - 1))), "[M-H]-")
      dt6 <- data.frame(cbind(
        "mz" = as.numeric(i.mz6), 
        "percent" = 20,
        "fragment" = paste0("[", input$sn1, "-H-H2O]-")))
      dt <- rbind(dt, dt6)
      
      i.mz7 <- MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)), 
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 - 
          (2*as.numeric(gsub(".*:", "", input$sn1))), 
        "N",
        "O", 2 - 1)))
      dt7 <- data.frame(cbind(
        "mz" = as.numeric(i.mz7), 
        "percent" = 5,
        "fragment" = paste0("[", input$sn1, "+N-O]-")))
      dt <- rbind(dt, dt7)
      
      i.mz8 <- MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)) - 2, 
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 - 
          (2*as.numeric(gsub(".*:", "", input$sn2))) - 5, 
        "O", 2 - 1)))
      dt8 <- data.frame(cbind(
        "mz" = as.numeric(i.mz8), 
        "percent" = 9,
        "fragment" = paste0("[", input$sn2, "-C2H5O]-")))
      dt <- rbind(dt, dt8)
      
      i.mz9 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)), 
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 - 
          (2*as.numeric(gsub(".*:", "", input$sn2))) - 2, 
        "O", 2 - 1))), "[M-H]-")
      dt9 <- data.frame(cbind(
        "mz" = as.numeric(i.mz9), 
        "percent" = 20,
        "fragment" = paste0("[", input$sn2, "-H-H2O]-")))
      dt <- rbind(dt, dt9)
      
      i.mz10 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C - 1, 
        "H", input$C*2 - (2*input$db) + 1 - 2, 
        "N",
        "O", 3 - 1))), "[M-H]-")
      dt10 <- data.frame(cbind("mz" = as.numeric(i.mz10), 
                               "percent" = 80,
                               "fragment" = "[M-H-CH2O]-"))
      dt <- rbind(dt, dt10)
    } else if(input$class == "CER_FA"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db) + 1, 
        "NO3"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M-H]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CER_hex"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6 - 6, 
        "H", input$C*2 - (2*input$db) + 11 - 10, 
        "N",
        "O", 3 + 5 - 5))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M-H-hexose]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CER_hex_FA"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, 
        "H", input$C*2 - (2*input$db) + 11, 
        "NO8"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M-H]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CER_lact"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6*2 - 6, 
        "H", input$C*2 - (2*input$db) + 1 + 10*2 - 10, 
        "N",
        "O", 3 + 5*2 - 5))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M-H-hexose]-"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6*2 - 6, 
        "H", input$C*2 - (2*input$db) + 1 + 10*2 - 10 - 2, 
        "N",
        "O", 3 + 5*2 - 5 - 1))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 75,
                              "fragment" = "[M-H-hexose-H2O]-"))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6*2 - 6*2, 
        "H", input$C*2 - (2*input$db) + 1 + 10*2 - 10*2, 
        "N",
        "O", 3 + 5*2 - 5*2))), "[M-H]-")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 80,
                              "fragment" = "[M-H-2(hexose)]-"))
      dt <- rbind(dt, dt3)
    } else if(input$class == "PA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)),
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 90,
                             "fragment" = paste0("[", input$sn2,"-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 5 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O6P"))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = paste0("[M-H-", input$sn1,"-H2O]-")))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 5 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1))) - 2), 
        "O7P"))), "[M-H]-")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 70,
                              "fragment" = paste0("[M-H-", input$sn1,"]-")))
      
      if(input$sn1 != input$sn2){
        i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn1)),
          "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
        dt4 <- data.frame(cbind("mz" = as.numeric(i.mz4), 
                                "percent" = 30,
                                "fragment" = paste0("[", input$sn1,"-H]-")))
        dt <- rbind(dt, dt4)
        
        i.mz5 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn2)), 
          "H", input$C*2 - (2*input$db) + 5 - 
            (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn2)))), 
          "O6P"))), "[M-H]-")
        dt5 <- data.frame(cbind("mz" = as.numeric(i.mz5), 
                                "percent" = 35,
                                "fragment" = paste0("[M-H-", input$sn2,"-H2O]-")))
        dt <- rbind(dt, dt5)
      } 
      dt <- rbind(dt, dt3)
    } else if(input$class == "mPA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)),
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn2, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 4 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2*input$db) + 7 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2), 
        "O7P"))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 20,
                              "fragment" = paste0("[M-H-", input$sn2, "]-")))
      
      if(input$sn1 != input$sn2){
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn1)),
          "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
        dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                                "percent" = 40,
                                "fragment" = paste0("[", input$sn1, "-H]-")))
        dt <- rbind(dt, dt3)
      } 
      dt <- rbind(dt, dt2)
    } else if(input$class == "dmPA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)),
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn1, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 5 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 9 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1))) - 2), 
        "O7P"))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 20,
                              "fragment" = paste0("[M-H-", input$sn1, "]-")))
      
      if(input$sn1 != input$sn2){
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn2)),
          "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
        dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                                "percent" = 35,
                                "fragment" = paste0("[", input$sn2, "-H]-")))
        dt <- rbind(dt, dt3)
      } 
      dt <- rbind(dt, dt2)
    } else if(input$class == "LPC"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7, "H", input$C*2 - (2*input$db) + 16, "NO7P"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M - CH3]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "LPC_MS3"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)),
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn1, "-H]-")))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "PC"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7, 
        "H", input$C*2 - ( 2*input$db) + 14, "NO8P"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M - CH3]-"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "LPE"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)),
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn1,"-H]-")))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "PE"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)),
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn2, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 5 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2*input$db) + 10 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2), 
        "NO7P"))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 20,
                              "fragment" = paste0("[M-H-", input$sn2, "]-")))
      
      if(input$sn1 != input$sn2){
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn1)),
          "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
        dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                                "percent" = 50,
                                "fragment" = paste0("[", input$sn1, "-H]-")))
        dt <- rbind(dt, dt3)
      } 
      dt <- rbind(dt, dt2)
    } else if(input$class == "PG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)),
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = paste0("[", input$sn2, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2)))), 
        "O8P"))), "[M-H]-")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 10,
        "fragment" = paste0("[M-H-", input$sn2, "-H2O]-")))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2), 
        "O9P"))), "[M-H]-")
      dt3 <- data.frame(cbind(
        "mz" = as.numeric(i.mz3), 
        "percent" = 15,
        "fragment" = paste0("[M-H-", input$sn2, "]-")))
      dt <- rbind(dt, dt3)
      
      i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn2)) - 3, 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2) - 8, 
        "O6P"))), "[M-H]-")
      dt4 <- data.frame(cbind(
        "mz" = as.numeric(i.mz4), 
        "percent" = 20,
        "fragment" = paste0("[M-H-", input$sn2, "-H2O-glycerol]-")))
      
      if(input$sn1 != input$sn2){
        i.mz5 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn1)),
          "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
        dt5 <- data.frame(cbind(
          "mz" = as.numeric(i.mz5), 
          "percent" = 70,
          "fragment" = paste0("[", input$sn1, "-H]-")))
        dt <- rbind(dt, dt5)
      } 
      i.mz6 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O8P"))), "[M-H]-")
      dt6 <- data.frame(cbind(
        "mz" = as.numeric(i.mz6), 
        "percent" = 25,
        "fragment" = paste0("[M-H-", input$sn1, "-H2O]-")))
      dt <- rbind(dt, dt6)
      
      i.mz7 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1))) - 2), 
        "O9P"))), "[M-H]-")
      dt7 <- data.frame(cbind(
        "mz" = as.numeric(i.mz7), 
        "percent" = 20,
        "fragment" = paste0("[M-H-", input$sn1, "]-")))
      dt <- rbind(dt, dt7)
      
      i.mz8 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - as.numeric(gsub(":.*", "", input$sn1)) - 3, 
        "H", input$C*2 - (2 + 2*input$db) + 13 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1))) - 2) - 8, 
        "O6P"))), "[M-H]-")
      dt8 <- data.frame(cbind(
        "mz" = as.numeric(i.mz8), 
        "percent" = 40,
        "fragment" = paste0("[M-H-", input$sn1, "-H2O-glycerol]-")))
      dt <- rbind(dt, dt8)
      dt <- rbind(dt, dt4)
    } else if(input$class == "LPI"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db), 
        "O2"))), "[M-H]-")
      dt <- data.frame(cbind(
        "mz" = as.numeric(i.mz1), 
        "percent" = 100,
        "fragment" = paste0("[", input$C, ":", input$db, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9 - 6, 
        "H", input$C*2 - (2*input$db) + 17 - 12, 
        "O", 12 - 6,
        "P"))), "[M-H]-")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 50,
        "fragment" = "[M-H-hexose-H2O]-"))
      
      dt <- rbind(dt, dt2)
    } else if(input$class == "PI"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn2)),
        "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 20,
                             "fragment" = paste0("[", input$sn2, "-H]-")))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2 + 2*input$db) + 17 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2)))), 
        "O11P"))), "[M-H]-")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 100,
        "fragment" = paste0("[M-H-", input$sn2, "-H2O]-")))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2 + 2*input$db) + 17 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2), 
        "O12P"))), "[M-H]-")
      dt3 <- data.frame(cbind(
        "mz" = as.numeric(i.mz3), 
        "percent" = 20,
        "fragment" = paste0("[M-H-", input$sn2, "]-")))
      dt <- rbind(dt, dt3)
      
      i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn2)) - 6, 
        "H", input$C*2 - (2 + 2*input$db) + 17 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2))) - 2) - 12, 
        "O6P"))), "[M-H]-")
      dt4 <- data.frame(cbind(
        "mz" = as.numeric(i.mz4), 
        "percent" = 70,
        "fragment" = paste0("[M-H-", input$sn2, "-H2O-inositol]-")))
      
      if(input$sn1 != input$sn2){
        i.mz5 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn1)),
          "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
        dt5 <- data.frame(cbind(
          "mz" = as.numeric(i.mz5), 
          "percent" = 50,
          "fragment" = paste0("[", input$sn1, "-H]-")))
        dt <- rbind(dt, dt5)
        
        i.mz6 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn1)), 
          "H", input$C*2 - (2 + 2*input$db) + 17 - 
            (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn1)))), 
          "O11P"))), "[M-H]-")
        dt6 <- data.frame(cbind(
          "mz" = as.numeric(i.mz6), 
          "percent" = 25,
          "fragment" = paste0("[M-H-", input$sn1, "-H2O]-")))
        dt <- rbind(dt, dt6)
        
        i.mz7 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn1)) - 6, 
          "H", input$C*2 - (2 + 2*input$db) + 17 - 
            (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn1))) - 2) - 12, 
          "O6P"))), "[M-H]-")
        dt7 <- data.frame(cbind(
          "mz" = as.numeric(i.mz7), 
          "percent" = 10,
          "fragment" = paste0("[M-H-", input$sn1, "-H2O-inositol]-")))
        dt <- rbind(dt, dt7)
        
        i.mz8 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn1)), 
          "H", input$C*2 - (2 + 2*input$db) + 17 - 
            (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn1))) - 2), 
          "O12P"))), "[M-H]-")
        dt8 <- data.frame(cbind(
          "mz" = as.numeric(i.mz8), 
          "percent" = 5,
          "fragment" = paste0("[M-H-", input$sn1, "]-")))
        dt <- rbind(dt, dt8)
      } 
      dt <- rbind(dt, dt4)
    } else if(input$class == "PS"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6 - 3, 
        "H", input$C*2 - (2*input$db) + 10 - 5, 
        "O8P"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M - H - serine]-"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 6 - 3 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 10 - 5 -
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O6P"))), "[M-H]-")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 20,
        "fragment" = paste0("[M - H - serine - ", input$sn1, "]-")))
      dt <- rbind(dt, dt2)
    } else if(input$class == "MGDG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 14, 
        "O10"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M-H]-"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", as.numeric(gsub(":.*", "", input$sn1)),
        "H", as.numeric(gsub(":.*", "", input$sn1))*2 -
          (2*(as.numeric(gsub(".*:", "", input$sn1)))), "O2"))), "[M-H]-")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = paste0("[", input$sn1,"-H]-")))
      dt <- rbind(dt, dt2)
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 14 -
          ((as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))) - 2), 
        "O9"))), "[M-H]-")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 20,
                              "fragment" = paste0("[M-H-", input$sn1,"]-")))
      if(input$sn1 == input$sn2){
        dt2$percent <- as.numeric(dt2$percent)*0.9
        dt3$percent <- as.numeric(dt3$percent)*0.3
      } else {
        dt2$percent <- as.numeric(dt2$percent)*0.5
        dt3$percent <- as.numeric(dt3$percent)*0.1
        
        i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", as.numeric(gsub(":.*", "", input$sn2)),
          "H", as.numeric(gsub(":.*", "", input$sn2))*2 -
            (2*(as.numeric(gsub(".*:", "", input$sn2)))), "O2"))), "[M-H]-")
        dt4 <- data.frame(cbind("mz" = as.numeric(i.mz4), 
                                "percent" = 20,
                                "fragment" = paste0("[", input$sn2,"-H]-")))
        dt <- rbind(dt, dt4)
        
        i.mz5 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn2)), 
          "H", input$C*2 - (2*input$db) + 14 -
            ((as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn2)))) - 2), 
          "O9"))), "[M-H]-")
        dt5 <- data.frame(cbind("mz" = as.numeric(i.mz5), 
                                "percent" = 15,
                                "fragment" = paste0("[M-H-", input$sn2,"]-")))
        dt <- rbind(dt, dt5)
      }
      
      dt <- rbind(dt, dt2)
      dt <- rbind(dt, dt3)
    } else if(input$class == "DGDG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 15, 
        "H", input$C*2 - (2*input$db) + 24, 
        "O15"))), "[M-H]-")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M-H]-"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 24 -
          ((as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1))))), 
        "O", 15 - 2))), "[M-H]-")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 100,
        "fragment" = paste0("[", input$sn1, "]-")))
      
      if(input$sn1 == input$sn2){
        dt2$percent <- as.numeric(dt2$percent)*0.4
      } else {
        dt2$percent <- as.numeric(dt2$percent)*0.2
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn2)), 
          "H", input$C*2 - (2*input$db) + 24 -
            ((as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn2))))), 
          "O", 15 - 2))), "[M-H]-")
        dt3 <- data.frame(cbind(
          "mz" = as.numeric(i.mz3), 
          "percent" = 20,
          "fragment" = paste0("[", input$sn2, "]-")))
        dt <- rbind(dt, dt3)
      }
      dt <- rbind(dt, dt2)
    }
    dt <- dt[order(dt$percent, decreasing = T),]
    dt <- dt[!duplicated(round(as.numeric(dt$mz), 4)),]
  })
  
  #### thr_MS2_POS ----
  thr_MS2_POS <- reactive({
    if(input$class == "CAR"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7 - 3, 
        "H", input$C*2 - (2 + 2*input$db) + 15 - 9, "O4"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M + H - trimethylamine]+"))
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db) - 2, "O"))), "[M+H]+")
      dt2 <- data.frame(cbind(
        "mz" = as.numeric(i.mz2), 
        "percent" = 20,
        "fragment" = paste0("[", input$C, ":", input$db, " + H - H2O]+")))
      dt <- rbind(dt, dt2)
    } else if(input$class == "CER"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, 
        "H", input$C*2 - (2*input$db) + 1 - 2, 
        "NO2"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M+H-H2O]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "CER_hex"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, 
        "H", input$C*2 - (2*input$db) + 11 - 2, 
        "N",
        "O", 3 + 5 - 1))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+H-H2O]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6 - 6, 
        "H", input$C*2 - (2*input$db) + 11 - 2 - 10, 
        "N",
        "O", 3 + 5 - 1 - 5))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 25,
                              "fragment" = "[M+H-hexose-H2O]+"))
      
      dt <- rbind(dt, dt2)
    } else if(input$class == "SM"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, 
        "H", input$C*2 - (2*input$db) + 13 - 2, 
        "N2O5P"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M+H-H2O]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "PA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, "H", input$C*2 - (2*input$db) + 5, 
        "O8P"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 30,
                             "fragment" = "[M+H]+"))
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 5 - 3, 
        "O4"))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = "[M+H-phosphate]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "mPA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 4, "H", input$C*2 - (2*input$db) + 7, 
        "O8P"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 30,
                             "fragment" = "[M+H]+"))
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 4 - 1, 
        "H", input$C*2 - (2*input$db) + 7 - 5, 
        "O4"))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = "[M+H-methylphosphate]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "dmPA"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 9, 
        "O8P"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 60,
                             "fragment" = "[M+H]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 5 - 2, 
        "H", input$C*2 - (2*input$db) + 9 - 7, 
        "O4"))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = "[M+H-dimethylphosphate]+"))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 5, 
        "H", input$C*2 - (2*input$db) + 9 - 3, 
        "O4"))), "[M+H]+")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 30,
                              "fragment" = "[M+H-phosphate]+"))
      dt <- rbind(dt, dt3)
    } else if(input$class == "LPC"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 18, 
        "NO6P"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+H-H2O]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula("C5H14NO4P")), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 30,
                              "fragment" = "[Phosphocholine]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "LPC_MS3"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8 - 3, 
        "H", input$C*2 - (2*input$db) + 18 - 9, 
        "O7P"))), "[M+Na]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M+Na-C3H9N]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "PC"){
      if(input$C < 36){
        i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
          "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 18, 
          "NO8P"))), "[M+H]+")
        dt <- data.frame(cbind("mz" = as.numeric(i.mz1 - 157.0504), 
                               "percent" = 100,
                               "fragment" = "[M + H - 157.0504]+"))
        
        i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn1)), 
          "H", input$C*2 - (2 + 2*input$db) + 18 - 
            (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn1)))), 
          "NO6P"
        ))), "[M+H]+")
        dt2 <- data.frame(cbind(
          "mz" = as.numeric(i.mz2), 
          "percent" = 40,
          "fragment" = paste0("[M + H - ", input$sn1, " - H2O]+")))
        
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn2)), 
          "H", input$C*2 - (2 + 2*input$db) + 18 - 
            (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn2)))) + 2, 
          "NO7P"))), "[M+H]+")
        dt3 <- data.frame(cbind(
          "mz" = as.numeric(i.mz3), 
          "percent" = 30,
          "fragment" = paste0("[M + H - ", input$sn2, "]+")))
        dt <- rbind(dt, dt3)
        if(input$C > 32){
          dt2$percent <- as.numeric(dt2$percent)*0.5
          
          i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn2)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn2)))), 
            "NO6P"
          ))), "[M+H]+")
          dt4 <- data.frame(cbind(
            "mz" = as.numeric(i.mz4), 
            "percent" = 10,
            "fragment" = paste0("[M + H - ", input$sn2, " - H2O]+")))
          dt <- rbind(dt, dt4)
        }
        dt <- rbind(dt, dt2)
      } else{
        if(input$sn1 == input$sn2){
          i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn1)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn1)))), 
            "NO6P"
          ))), "[M+H]+")
          dt2 <- data.frame(cbind(
            "mz" = as.numeric(i.mz2), 
            "percent" = 100,
            "fragment" = paste0("[M + H - ", input$sn1, " - H2O]+")))
          
          i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn2)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn2)))) + 2, 
            "NO7P"))), "[M+H]+")
          dt3 <- data.frame(cbind(
            "mz" = as.numeric(i.mz3), 
            "percent" = 80,
            "fragment" = paste0("[M + H - ", input$sn2, "]+")))
          dt <- rbind(dt2, dt3)
        } else {
          i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn2)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn2)))) + 2, 
            "NO7P"))), "[M+H]+")
          dt3 <- data.frame(cbind(
            "mz" = as.numeric(i.mz3), 
            "percent" = 100,
            "fragment" = paste0("[M + H - ", input$sn2, "]+")))
          
          i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn1)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn1)))), 
            "NO6P"
          ))), "[M+H]+")
          dt2 <- data.frame(cbind(
            "mz" = as.numeric(i.mz2), 
            "percent" = 80,
            "fragment" = paste0("[M + H - ", input$sn1, " - H2O]+")))
          dt <- rbind(dt2, dt3)
          
          i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
            "C", input$C + 8 - as.numeric(gsub(":.*", "", input$sn2)), 
            "H", input$C*2 - (2 + 2*input$db) + 18 - 
              (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
                gsub(".*:", "", input$sn2)))), 
            "NO6P"
          ))), "[M+H]+")
          dt4 <- data.frame(cbind(
            "mz" = as.numeric(i.mz4), 
            "percent" = 40,
            "fragment" = paste0("[M + H - ", input$sn2, " - H2O]+")))
          dt <- rbind(dt, dt4)
        }
      }
    } else if(input$class == "LPE"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, "H", input$C*2 - (2 + 2*input$db) + 6, 
        "O3"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 30,
                             "fragment" = "[M+H-phosphoethanolamina]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2 + 2*input$db) + 12, "NO6P"))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 100,
                              "fragment" = "[M+H-H2O]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "PE"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, "H", input$C*2 - (2*input$db) + 2, 
        "O4"))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M+H-phosphoethanolamine]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "PG"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, 
        "H", input$C*2 - (2*input$db) + 11, "O10P"))), 
        "[M+NH4]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz - 141.04), 
                             "percent" = 100,
                             "fragment" = "[M+NH4-C3H9O6]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "LPI"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 17 - 2, 
        "O", 12 - 1, 
        "P"))), "[M+H]+") 
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+H-H2O]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 17 - 2*2, 
        "O", 12 - 1*2, 
        "P"))), "[M+H]+") 
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 30,
                              "fragment" = "[M+H-2(H2O)]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "PI"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2 + 2*input$db) + 17, "O13P"))), "[M+NH4]+") - 
        MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) + 0.984
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+NH4-phosphoinositol]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 15, 
        "O13P"))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 60,
                              "fragment" = "[M+H]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "PS"){
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6 - 3, 
        "H", input$C*2 - (2*input$db) + 10 - 8, 
        "O", 10 - 6))), 
        "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz), 
                             "percent" = 100,
                             "fragment" = "[M + H - phosphoserine]+"))
      dt <- rbind(dt, dt)
      dt <- dt[1,]
    } else if(input$class == "MGDG"){
      i.mz1 <- as.numeric(mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - 6, 
        "H", input$C*2 - (2*input$db) + 14 - 12, 
        "O", 10 - 6))), "[M+NH4]+")) + 0.984
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+NH4-galactose]+ +0.984"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 9 - 6, 
        "H", input$C*2 - (2*input$db) + 14 - 12, 
        "O", 10 - 6))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 20,
                              "fragment" = "[M+H-galactose]+"))
      dt <- rbind(dt, dt2)
      
    } else if(input$class == "MGDG_Na"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 14 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O", 10 - 2))), "[M+Na]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = paste0("[M+Na-", input$sn1, "]+")))
      #dt <- rbind(dt, dt)
      #dt <- dt[1,]
      
      #fml2 <- "H2O"
      #dt2 <- IsotopicDistribution(formula = ListFormula(fml2))
      #dt2[,3] <- dt2[,3]*0.05
      i.mz2 <- mass2mz(MonoisotopicMass(formula = ListFormula("H2O")), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 5,
                              "fragment" = "X"))
      #dt <- rbind(dt, dt2)
      
      if(input$sn1 != input$sn2){
        i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
          "C", input$C + 9 - as.numeric(gsub(":.*", "", input$sn2)), 
          "H", input$C*2 - (2*input$db) + 14 - 
            (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
              gsub(".*:", "", input$sn2)))), 
          "O", 10 - 2))), "[M+Na]+")
        dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                                "percent" = 50,
                                "fragment" = paste0("[M+Na-", input$sn2, "]+")))
        dt <- rbind(dt, dt3)
      }
      dt <- rbind(dt, dt2)
      dt <- dt[dt$mz > 20, ]
    } else if(input$class == "DGDG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - 12, 
        "H", input$C*2 - (2*input$db) + 24 - 22, 
        "O", 15 - 11))), "[M+NH4]+") + 0.984
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+NH4-2(galactose)]+ +0.984"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - 12, 
        "H", input$C*2 - (2*input$db) + 24 - 22, 
        "O", 15 - 11))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 20,
                              "fragment" = "[M+H-2(galactose)]+"))
      dt <- rbind(dt, dt2)
    } else if(input$class == "DGDG_Na"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2*input$db) + 24 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2)))), 
        "O", 15 - 2))), "[M+Na]+")
      dt <- data.frame(cbind(
        "mz" = as.numeric(i.mz1), 
        "percent" = 100,
        "fragment" = paste0("[M+Na-", input$sn2, "]+")))
      
      #i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
      #  "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn1)) - 6, 
      #  "H", input$C*2 - (2*input$db) + 24 - 
      #    (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
      #      gsub(".*:", "", input$sn1)))) - 10, 
      #  "O", 15 - 2 - 5))), "[M+Na]+")
      #dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
      #                        "percent" = 20,
      #                        "fragment" = "X"))
      #dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 24 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O", 15 - 2))), "[M+Na]+")
      dt3 <- data.frame(cbind(
        "mz" = as.numeric(i.mz3), 
        "percent" = 60,
        "fragment" = paste0("[M+Na-", input$sn1, "]+")))
      dt <- rbind(dt, dt3)
      
      i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - 6, 
        "H", input$C*2 - (2*input$db) + 24 - 10, 
        "O", 15 - 5))), "[M+Na]+")
      dt4 <- data.frame(cbind(
        "mz" = as.numeric(i.mz4), 
        "percent" = 20,
        "fragment" = "[M+Na-hexose]+"))
      dt <- rbind(dt, dt4)
      
      i.mz5 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 15 - as.numeric(gsub(":.*", "", input$sn2)) - 6, 
        "H", input$C*2 - (2*input$db) + 24 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2)))) - 10, 
        "O", 15 - 2 - 5))), "[M+Na]+")
      dt5 <- data.frame(cbind(
        "mz" = as.numeric(i.mz5), 
        "percent" = 15,
        "fragment" = paste0("[M+Na-", input$sn2, "-hexose]+")))
      dt <- rbind(dt, dt5)
      
    } else if(input$class == "DAG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
        paste0("C", input$C + 3, 
               "H", input$C*2 - (2*input$db) + 4, 
               "O", 5))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+H]+"))
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(
        paste0("C", input$C + 3, 
               "H", input$C*2 - (2*input$db) + 4 - 2, 
               "O", 5 - 1))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 30,
                              "fragment" = "[M+H-H2O]+"))
      dt <- rbind(dt, dt2)
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 4 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O", 5 - 2))), "[M+H]+")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 20,
                              "fragment" = paste0("[M+H-", input$sn1,"]+")))
      dt <- rbind(dt, dt3)
    } else if(input$class == "TAG"){
      i.mz1 <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 2, 
        "O", 6))), "[M+H]+")
      dt <- data.frame(cbind("mz" = as.numeric(i.mz1), 
                             "percent" = 100,
                             "fragment" = "[M+H]+"))
      
      i.mz2 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn1)), 
        "H", input$C*2 - (2*input$db) + 2 - 
          (as.numeric(gsub(":.*", "", input$sn1))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn1)))), 
        "O", 6 - 2))), "[M+H]+")
      dt2 <- data.frame(cbind("mz" = as.numeric(i.mz2), 
                              "percent" = 50,
                              "fragment" = paste0("[M+H-", input$sn1, "]+")))
      dt <- rbind(dt, dt2)
      
      i.mz3 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn2)), 
        "H", input$C*2 - (2*input$db) + 2 - 
          (as.numeric(gsub(":.*", "", input$sn2))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn2)))), 
        "O", 6 - 2))), "[M+H]+")
      dt3 <- data.frame(cbind("mz" = as.numeric(i.mz3), 
                              "percent" = 70,
                              "fragment" = paste0("[M+H-", input$sn2, "]+")))
      dt <- rbind(dt, dt3)
      
      i.mz4 <- mass2mz(MonoisotopicMass(formula =ListFormula(paste0(
        "C", input$C + 3 - as.numeric(gsub(":.*", "", input$sn3)), 
        "H", input$C*2 - (2*input$db) + 2 - 
          (as.numeric(gsub(":.*", "", input$sn3))*2 - (2*as.numeric(
            gsub(".*:", "", input$sn3)))), 
        "O", 6 - 2))), "[M+H]+")
      dt4 <- data.frame(cbind("mz" = as.numeric(i.mz4), 
                              "percent" = 30,
                              "fragment" = paste0("[M+H-", input$sn3, "]+")))
      dt <- rbind(dt, dt4)
    } 
    dt <- dt[order(dt$percent, decreasing = T),]
    dt <- dt[!duplicated(round(as.numeric(dt$mz), 4)),]
  })
  
  #### thr_MS2_N ----
  output$thr_MS2_N <- renderPlot({
    dt <- thr_MS2_NEG()
    dt$mz <- as.numeric(dt$mz)
    dt$percent <- as.numeric(dt$percent)
    if(input$class == "FFA"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db), "O2"))), ad)
    } else if(input$class == "FFA_red"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db) + 2, "O"))), ad)
    } else if(input$class == "FFA_PA"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, 
        "H", input$C*2 - (2*input$db) + 8, "O", 2 + 2))), ad)
    } else if(input$class == "CAR"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7, "H", input$C*2 - (2 + 2*input$db) + 15, "NO4"))), ad)
    } else if(input$class == "CER"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db) + 1, 
        "NO3"))), ad)
    } else if(input$class == "CER_FA"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db) + 1, 
        "NO3"))), ad)
    } else if(input$class == "CER_hex"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, "H", input$C*2 - (2*input$db) + 11, 
        "NO8"))), ad)
    } else if(input$class == "CER_hex_FA"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, "H", input$C*2 - (2*input$db) + 11, 
        "NO8"))), ad)
    } else if(input$class == "CER_lact"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6*2, "H", input$C*2 - (2*input$db) + 1 + 10*2, 
        "N", "O", 3 + 5*2))), ad)
    } else if(input$class == "SM"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 13, 
        "N2O6P"))), ad)
    } else if(input$class == "PA"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, "H", input$C*2 - (2*input$db) + 5, "O8P"))), ad)
    } else if(input$class == "mPA"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 4, "H", input$C*2 - (2*input$db) + 7, "O8P"))), ad)
    } else if(input$class == "dmPA"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 9, "O8P"))), ad)
    } else if(input$class == "LPC"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 20, "NO7P"))), ad)
    } else if(input$class == "LPC_MS3"){
      ad <- "[M+CHO2]-"
      i.mz <- MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7, "H", input$C*2 - (2 + 2*input$db) + 17, "NO7P")))
    } else if(input$class == "PC"){
      ad <- "[M+CHO2]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 18, 
        "NO8P"))), ad)
    } else if(input$class == "LPE"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2 + 2*input$db) + 14, "NO7P"))), ad)
    }else if(input$class == "PE"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 10, "NO8P"))), ad)
    } else if(input$class == "PG"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, "H", input$C*2 - (2 + 2*input$db) + 13, "O10P"))), ad)
    } else if(input$class == "LPI"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 17, "O12P"))), ad)
    } else if(input$class == "PI"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 17, "O13P"))), ad)
    } else if(input$class == "PS"){
      ad <- "[M-H]-"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, 
        "H", input$C*2 - (2*input$db) + 10, "NO10P"))), ad)
    } else if(input$class == "MGDG" | input$class == "MGDG_Na"){
      ad <- "[M+CHO2]-"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 14, "O10"))), ad)
    } else if(input$class == "DGDG" | input$class == "DGDG_Na"){
      ad <- "[M+CHO2]-"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 15, 
        "H", input$C*2 - (2*input$db) + 24, "O15"))), ad)
    } else if(input$class == "DAG"){
      ad <- "[M+CHO2]-"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 4, "O5"))), ad)
    } else if(input$class == "TAG"){
      ad <- "[M+CHO2]-"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 2, "O6"))), ad)
    } 
    plot(dt$mz, dt$percent,
         type = "h", xlim = c(50, i.mz), 
         ylim = c(0, 105),
         xlab = "m/z", ylab = "relative intensity", 
         main = paste("MS2", ad, sprintf("%.5f", i.mz)))
    idx <- which(dt$percent > 15)
    text(dt$mz[idx],
         dt$percent[idx], 
         sprintf("%.5f", dt$mz[idx]), pos = 3)
  })
  
  #### thr_MS2_P ----
  output$thr_MS2_P <- renderPlot({
    dt <- thr_MS2_POS()
    dt$mz <- as.numeric(dt$mz)
    dt$percent <- as.numeric(dt$percent)
    if(input$class == "CAR"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 7, "H", input$C*2 - (2 + 2*input$db) + 15, "NO4"))), ad)
    } else if(input$class == "CER"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C, "H", input$C*2 - (2*input$db) + 1, 
        "NO3"))), ad)
    } else if(input$class == "CER_hex"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, "H", input$C*2 - (2*input$db) + 11, 
        "NO8"))), ad)
    } else if(input$class == "SM"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 13, 
        "N2O6P"))), ad)
    } else if(input$class == "PA"){
      ad <- "[M+NH4]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, "H", input$C*2 - (2*input$db) + 5, "O8P"))), ad)
    } else if(input$class == "mPA"){
      ad <- "[M+NH4]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 4, "H", input$C*2 - (2*input$db) + 7, "O8P"))), ad)
    } else if(input$class == "dmPA"){
      ad <- "[M+NH4]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 9, "O8P"))), ad)
    } else if(input$class == "LPC"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 20, "NO7P"))), ad)
    } else if(input$class == "LPC_MS3"){
      ad <- "[M+Na]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 20, "NO7P"))), ad)
    } else if(input$class == "PC"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 8, "H", input$C*2 - (2 + 2*input$db) + 18, 
        "NO8P"))), ad)
    } else if(input$class == "LPE"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2 + 2*input$db) + 14, "NO7P"))), ad)
    }else if(input$class == "PE"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 5, "H", input$C*2 - (2*input$db) + 10, "NO8P"))), ad)
    } else if(input$class == "PG"){
      ad <- "[M+NH4]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, "H", input$C*2 - (2 + 2*input$db) + 13, "O10P"))), ad)
    } else if(input$class == "LPI"){
      ad <- "[M+H]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 17, "O12P"))), ad)
    } else if(input$class == "PI"){
      ad <- "[M+NH4]+"
      i.mz <- mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2 + 2*input$db) + 17, "O13P"))), ad)
    } else if(input$class == "PS"){
      ad <- "[M+H]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 6, 
        "H", input$C*2 - (2*input$db) + 10, "NO10P"))), ad)
    } else if(input$class == "MGDG"){
      ad <- "[M+NH4]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 14, "O10"))), ad)
    } else if(input$class == "MGDG_Na"){
      ad <- "[M+Na]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 9, 
        "H", input$C*2 - (2*input$db) + 14, "O10"))), ad)
    } else if(input$class == "DGDG"){
      ad <- "[M+NH4]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 15, 
        "H", input$C*2 - (2*input$db) + 24, "O15"))), ad)
    } else if(input$class == "DGDG_Na"){
      ad <- "[M+Na]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 15, 
        "H", input$C*2 - (2*input$db) + 24, "O15"))), ad)
    } else if(input$class == "DAG"){
      ad <- "[M+NH4]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 4, "O5"))), ad)
    } else if(input$class == "TAG"){
      ad <- "[M+NH4]+"
      i.mz <-  mass2mz(MonoisotopicMass(formula = ListFormula(paste0(
        "C", input$C + 3, 
        "H", input$C*2 - (2*input$db) + 2, "O6"))), ad)
    } 
    plot(dt$mz, dt$percent,
         type = "h", xlim = c(50, i.mz), 
         ylim = c(0, 105),
         xlab = "m/z", ylab = "relative intensity", 
         main = paste("MS2", ad, sprintf("%.5f", i.mz)))
    idx <- which(dt$percent > 15)
    text(dt$mz[idx],
         dt$percent[idx], 
         sprintf("%.5f", dt$mz[idx]), pos = 3)
  })
  
  output$thr_MS2_tb_N <- renderPrint({
    dt <- thr_MS2_NEG()
    dt$mz <- sprintf("%.6f", as.numeric(dt$mz))
    dt <- dt[order(dt$mz),]
    dt
  })
  
  output$thr_MS2_tb_P <- renderPrint({
    dt <- thr_MS2_POS()
    dt$mz <- sprintf("%.6f", as.numeric(dt$mz))
    dt <- dt[order(dt$mz),]
    dt
  })
  
} # close server

shinyApp(ui = ui, server = server)
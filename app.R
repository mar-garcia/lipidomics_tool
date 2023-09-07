options(repos = BiocManager::repositories())
options("repos")

cmps <- readxl::read_xlsx("compounds.xlsx")

library(shiny)
library(MetaboCoreUtils)
library(plotly)
library(tidyverse)


.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

# to add mPA, dmPA
source("lipid_functions.R")

load("lipid_workspace.RData")

cmps_db$pos <- cmps_db$neg <- NA
idx <- which(cmps_db$class %in% c("FA", "CAR", 
                                  "SM", "Cer", "Cer;O3", "Cer;O4", 
                                  "HexCer", "HexCer;O3", "HexCer;O4", "LactCer", 
                                  "LPC", "LPE", "LPS", "PC", "PE", "PS",
                                  "DGTS"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+H]+")
idx <- which(cmps_db$class %in% c("LPA", "LPG", "LPI",
                                  "PA", "PG", "PI", "MGDG", "DGDG", "DGGA", 
                                  "MG", "DG", "TG"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+NH4]+")
idx <- which(cmps_db$class %in% c("pHexFA"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+Na]+")
idx <- which(cmps_db$class %in% c("FA", "LPA", "LPE", "LPG", "LPI", "LPS", 
                                  "PA", "PE", "PG", "PI", "PS", 
                                  "MG", "DG", "TG", "DGGA"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M-H]-")
idx <- which(cmps_db$class %in% c("pHexFA", "CAR", 
                                  "SM", "Cer", "Cer;O3", "Cer;O4", 
                                  "HexCer", "HexCer;O3", "HexCer;O4", "LactCer", 
                                  "LPC", "PC", "MGDG", "DGDG", "DGTS"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M+CHO2]-")



mz_calculator <- function(class, fml){
  if(class == "FA"){
    mass2mz(calculateMass(fml), c("[M-H]-"))
  } else if(class %in% c("MG", "DG", "TG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+"))
  } else if(class %in% c("LPE", "LPS", "PE", "PS")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M-H]-"))
  } else if(class %in% c("Cer", "Cer;O3", "Cer;O4", 
                         "HexCer", "HexCer;O3", "HexCer;O4", "LactCer",
                         "SM", "LPC", "PC", "DGTS")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M+CHO2]-"))
  } else if(class %in% c("CAR", "LPA", "LPG", "LPI", "PA", "PG", "PI", "DGGA")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M-H]-"))
  } else if(class %in% c("MGDG", "DGDG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M+CHO2]-"))
  } else if(class %in% c("pHexFA")){
    mass2mz(calculateMass(fml), c("[M+Na]+", "[M+CHO2]-"))
  }
}

ref <- c(12, 1, 16, 14, 31, 32)
names(ref) <- c("C", "H", "O", "N", "P", "S")

cmps$formula <- fml_maker(cmps$class, cmps$C, cmps$db)
cmps$mass <- calculateMass(cmps$formula)

cmps$NM <- NA
for(i in seq(nrow(cmps))){
  myform <- countElements(cmps$formula[i])
  pippo <- ref[names(myform[[1]])]
  cmps$NM[i] <- t(pippo) %*% myform[[1]]
}

cmps$KMD <- cmps$NM - (cmps$mass * ref["H"] / calculateMass("H"))
cmps$db <- as.character(cmps$db)


sn <- data.frame(
  "C" = rep(seq(3, 25), each = 6),
  "db" = rep(seq(0, 5), (25-3+1))
)
sn$formula <- paste0("C", sn$C, "H", sn$C*2 - 2*sn$db, "O2")
sn <- sn[!grepl("-", sn$formula),]
sn$sn <- paste0(sn$C, ":", sn$db)
sn$mass <- NA
sn.list <- list()
for(i in seq(nrow(sn))){
  sn$mass[i] <- calculateMass(sn$formula[i])
  sn.list[[i]] <- sn$sn[i]
  names(sn.list)[[i]] <- sn$sn[i]
}

mzdif.pos1 <- data.frame(rbind(
  c(calculateMass("NH3"), "loss NH3 -> LPA / PA / PG / MGDG / MG / DG / TG"),
  c(calculateMass("NH3H2O"), "loss NH3 & H2O -> LPA / MGDG / DG / DGGA"),
  c(calculateMass("H3PO4NH3"), "loss NH3 & PA -> PA"),
  c(calculateMass("CH3"), "loss CH3"),
  cbind(sn$mass + calculateMass("NH3"), paste("loss NH3 &", sn$sn, "-> DG / TG")),
  c(calculateMass("H2O"), "loss H2O -> LPC / LPS"),
  c(calculateMass("C6H10O5"), "loss hexose -> HexCer"),
  c(calculateMass("C6H10O5H2O"), "loss hexose & H2O -> HexCer"),
  c(calculateMass("C6H10O5H2OH2O"), "loss hexose & 2(H2O) -> HexCer"),
  c(mass2mz(calculateMass("C5H14NO4P"), "[M+H]+"), "[PC+H]+ -> LPC"),
  c(calculateMass("C3H7NO3"), "loss serine -> LPS"),
  c(calculateMass("C3H7NO3HPO3"), "loss PS -> LPS / PS"),
  c(calculateMass("C2H8NO4P"), "loss PE -> PE"),
  c(calculateMass("C3H8O3HPO3NH3"), "loss NH3 & PG -> PG"),
  c(calculateMass("C6H10O5NH3"), "loss NH3 & hexose -> MGDG"),
  c(calculateMass("C6H12O6NH3"), "loss NH3 & hexose & H2O -> MGDG"),
  c(calculateMass("C6H10O5C6H10O5NH3"), "loss NH3 & 2(hexose) -> DGDG"),
  c(calculateMass("C6H10O5C6H12O6NH3"), "loss NH3 & 2(hexose) & H2O -> DGDG"),
  c(mass2mz(calculateMass(subtractElements(addElements(addElements("C4H9NO3", "CH3CH3CH3"), "C3H8O3"), "H2OH3"))), "[GTS+H]+ -> DGTS"),
  c(calculateMass("NH3C6H10O7"), "loss NH3 & glucuronide & H2O -> DGGA"),
  c(calculateMass("NH3C6H8O6"), "loss NH3 & glucuronide -> DGGA"),
  
  cbind(sn$mass, paste0("loss '", sn$sn, "-H2O' -> PC / DGTS / TG;O / TG;O2")),
  cbind(sn$mass - calculateMass("H2O"), paste("loss ", sn$sn, " -> PC/ DGTS"))
  
))
mzdif.pos2 <- data.frame(rbind(
  cbind(mass2mz(calculateMass(addElements(sn$formula, "C3H4O")), "[M+H]+"), paste0("[", sn$sn, "+H+C3H4O]+ -> MGDG / DGDG / DGGA")),
  cbind(calculateMass(subtractElements(addElements(sn$formula, "NH2"), "O")), paste0("[", sn$sn, "+NH2-O]+ -> HexCer")),
  cbind(calculateMass(subtractElements(addElements(sn$formula, "N"), "O2")), paste0("[", sn$sn, "+N-O2]+ -> HexCer")),
  cbind(calculateMass(addElements(sn$formula, "C3H4N")), paste0("[", sn$sn, "+C3H4N]+ -> HexCer")),
  cbind(calculateMass(addElements(sn$formula, "C3H6NO")), paste0("[", sn$sn, "+C3H6NO]+ -> HexCer")),
  cbind(mass2mz(calculateMass(addElements(sn$formula, "O")), "[M+Na]+"), paste0("[", sn$sn, ";O+Na]+ -> TG;O")),
  cbind(mass2mz(calculateMass(addElements(sn$formula, "O2")), "[M+Na]+"), paste0("[", sn$sn, ";O2+Na]+ -> TG;O2"))
  
))
colnames(mzdif.pos1) <- colnames(mzdif.pos2) <- c("dif", "add")
mzdif.pos1$dif <- as.numeric(mzdif.pos1$dif)
mzdif.pos2$dif <- as.numeric(mzdif.pos2$dif)

mzdif.neg1 <- data.frame(rbind(
  c(calculateMass("HCOOH"), "loss HCOOH -> pHexFA / Cer / HexCer / MGDG / DGDG"),
  c(calculateMass("HCOOHCH2"), "loss HCOOH & CH3 -> SM / LPC / PC"),
  c(calculateMass("H2O"), "loss H2O -> Cer / H2O"),
  c(calculateMass("CH2O"), "loss CH2O -> Cer"),
  c(calculateMass("CH2OH2O"), "loss CH2O & H2O -> Cer"),
  c(calculateMass("CH3OH"), "loss CH3OH -> Cer"),
  c(calculateMass("C6H10O5"), "loss hexose -> HexCer;O3"),
  c(calculateMass("C6H12O6"), "loss hexose & H2O -> HexCer;O3"),
  c(calculateMass(subtractElements("C3H7NO3", "H2O")), "loss Ser -> LPS / PS"),
  c(calculateMass("C4H11NO"), "loss C4H11NO -> SM"),
  
  cbind(sn$mass, paste("loss", sn$sn, "& H2O -> LPI / PA / PI / DGGA")),
  cbind(sn$mass - calculateMass("H2O"), paste("loss", sn$sn, "-> PA / PE / PG / PI / DGGA")),
  cbind(sn$mass - calculateMass("H2O") + calculateMass("HCOOH"), paste("loss HCOOH &", sn$sn, "-> MGDG")),
  cbind(sn$mass - calculateMass("H2O") + calculateMass("C3H8O3"), paste("loss", sn$sn, "& glycerol -> PG")),
  cbind(sn$mass - calculateMass("H2O") + calculateMass("C6H12O6"), paste("loss", sn$sn, "& inositol -> PI")),
  cbind(calculateMass(addElements("C3H5NO2", sn$formula)), paste0("loss Ser & ", sn$sn, " -> PS"))
))
mzdif.neg2 <- data.frame(rbind(
  cbind(mass2mz(sn$mass, "[M-H]-"), paste0("[", sn$sn, "-H]- -> pHexFA / Cer (sn1) / LPG / LPI / PA / PE / PG / PI / MGDG / DGGA")),
  cbind(mass2mz(calculateMass(addElements(sn$formula, "C2H3N")), "[M-H]-"), paste0("[", sn$sn, "-H+C2H3N]- -> Cer (Ssn1) / HexCer;O3 / LPI")),
  cbind(mass2mz(calculateMass(subtractElements(addElements(sn$formula, "NH3C2"), "O")), "[M-H]-"), paste0("[", sn$sn, "-OH+C2H3N]- -> Cer (Tsn1)")),
  cbind(mass2mz(calculateMass(subtractElements(addElements(sn$formula, "NH3C2"), "CO")), "[M-H]-"), paste0("[", sn$sn, "-CHO+C2H3N]- -> Cer (Xsn1)")),
  #cbind(mass2mz(calculateMass(subtractElements(sn$formula, "H2O")), "[M-H]-"), paste0("[", sn$sn, "-H-H2O]- -> Cer (Rsn2)")),
  cbind(mass2mz(calculateMass(addElements(subtractElements(sn$formula, "O"), "NH")), "[M-H]-"), paste0("[", sn$sn, "-OH+NH]- -> Cer (Usn1)")),
  #cbind(mass2mz(calculateMass(subtractElements(sn$formula, "C2H4O")), "[M-H]-"), paste0("[", sn$sn, "-H-C2H4O]- -> Cer (Psn2)")),
  #cbind(calculateMass(subtractElements(addElements(sn$formula, "N"), "H2OCH2")), paste0("[M-H-", sn$sn, "(sn-N)] -> SM")),
  #cbind(calculateMass(subtractElements(sn$formula, "H2OCH2")) + calculateMass("C5H13NO"), paste0("[M-H-C5H13NO-", sn$sn, "(sn-N)] -> SM")),
  
  cbind(mass2mz(calculateMass(subtractElements("C3H9O6P", "H2O")), "[M-H]-"), "[GP-H-H2O]- -> LPA"),
  cbind(mass2mz(calculateMass(subtractElements("C6H13O9P", "H2O")), "[M-H]-"), "[IP-H-H2O]- -> LPI")
))

colnames(mzdif.neg1) <- colnames(mzdif.neg2) <- c("dif", "add")
mzdif.neg1$dif <- as.numeric(mzdif.neg1$dif)
mzdif.neg2$dif <- as.numeric(mzdif.neg2$dif)

# UI ------------------------------------------------------------------------
ui <- navbarPage(
  "Lipidomics",
  
  tabPanel(
    "Main Panel",
    column(6, h3("Formula"),
           fluidRow(
             column(2, selectInput(inputId = "class", label = "Lipid class:",
                                   choices = list("FA" = "FA",
                                                  "pHexFA" = "pHexFA",
                                                  "CAR" = "CAR", 
                                                  "SM" = "SM",
                                                  "Cer" = "Cer",
                                                  "Cer;O3" = "Cer;O3",
                                                  "Cer;O4" = "Cer;O4",
                                                  "HexCer" = "HexCer",
                                                  "HexCer;O3" = "HexCer;O3",
                                                  "HexCer;O4" = "HexCer;O4",
                                                  "LactCer" = "LactCer",
                                                  "LPA" = "LPA",
                                                  "LPC" = "LPC",
                                                  "LPE" = "LPE",
                                                  "LPG" = "LPG",
                                                  "LPI" = "LPI",
                                                  "LPS" = "LPS",
                                                  "PA" = "PA",
                                                  "PC" = "PC",
                                                  "PE" = "PE",
                                                  "PG" = "PG",
                                                  "PI" = "PI",
                                                  "PS" = "PS",
                                                  "MGDG" = "MGDG",
                                                  "DGDG" = "DGDG",
                                                  "DGGA" = "DGGA",
                                                  "DGTS" = "DGTS",
                                                  "MG" = "MG",
                                                  "DG" = "DG",
                                                  "TG" = "TG"))),
             column(2, numericInput(inputId = "C", label = "C", value = 18)),
             column(2, numericInput(inputId = "db", label = "db", value = 0))
           ),
           column(4, verbatimTextOutput("formula")),
           fluidRow(),
           fluidRow(
             h3("Major adducts"),
             column(6, fluidRow(verbatimTextOutput("mzvals"))))
    ), # close left column
    column(6, h3("Compound Prediction"),
           column(3, 
                  numericInput("mzs", "m/z", value = 283.2643),
                  radioButtons("pol", label = "Polarity:",
                               choices = list("ESI+" = "pos", 
                                              "ESI-" = "neg"),
                               selected = "neg"),
                  numericInput("ppms", "ppm", value = 5)
           ),
           column(3, verbatimTextOutput("cmps_pred"))
    ), # close right column
    fluidRow(),
    h3("Theoretical MS2"),
    column(4, h2("ESI+"),
           fluidRow(plotOutput("ms2_thr_pos"))),
    column(1),
    column(4, h2("ESI-"),
           fluidRow(plotOutput("ms2_thr_neg"))),
    column(1, selectInput("sn1", "sn1",
                          choices = sn$sn,
                          selected = "18:0")),
    column(1, selectInput("sn2", "sn2",
                          choices = sn$sn,
                          selected = "18:0")),
    column(1, selectInput("sn3", "sn3",
                          choices = sn$sn,
                          selected = "18:0"))
  ), # close tabPanel "P1"
  #tabPanel(
  #  "Kendrick Mass Defect (KMD)",
  #  sidebarLayout(
  #    sidebarPanel(
  #      checkboxGroupInput("kmd_class", label = "Lipid classes:", 
  #                         choices = list("FA" = "FA",
  #                                        "pHexFA" = "pHexFA",
  #                                        "CAR" = "CAR", 
  #                                        "SM" = "SM",
  #                                        "Cer" = "Cer",
  #                                        "Cer;O3" = "Cer;O3",
  #                                        "Cer;O4" = "Cer;O4",
  #                                        "HexCer" = "HexCer",
  #                                        "HexCer;O3" = "HexCer;O3",
  #                                        "HexCer;O4" = "HexCer;O4",
  #                                        "LactCer" = "LactCer",
  #                                        "LPA" = "LPA",
  #                                        "LPC" = "LPC",
  #                                        "LPE" = "LPE",
  #                                        "LPG" = "LPG",
  #                                        "LPI" = "LPI",
  #                                        "LPS" = "LPS",
  #                                        "PA" = "PA",
  #                                        "PC" = "PC",
  #                                        "PE" = "PE",
  #                                        "PG" = "PG",
  #                                        "PI" = "PI",
  #                                        "PS" = "PS",
  #                                        "MGDG" = "MGDG",
  #                                        "DGDG" = "DGDG",
  #                                        "DGTS" = "DGTS",
  #                                        "MG" = "MG",
  #                                        "DG" = "DG",
  #                                        "TG" = "TG"),
  #                         selected = "TG"),
  #      radioButtons("kmd_color", label = "Color by:",
  #                   choices = list("Class" = "class", 
  #                                  "Double bonds" = "db"),
  #                   selected = "class")
  #    ),
  #    mainPanel(
  #      plotlyOutput("kmd")
  #    )
  #  )), # close tabPanel "KMD"
  tabPanel(
    "MS2",
    column(6, h2("ESI+"),
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
           ),
           fluidRow(
             column(2, numericInput("posfrag5", "Fragment 5", value = 0)),
             column(8, verbatimTextOutput("posfrag5add"))
           ),
           fluidRow(
             column(2, numericInput("posfrag6", "Fragment 6", value = 0)),
             column(8, verbatimTextOutput("posfrag6add"))
           )), # close column ESI+
    column(6, h2("ESI-"),
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
           )) # close column ESI-
  ) # close tabPanel MS2
)

# SERVER ---------------------------------------------------------------------
server <- function(input, output) {
  
  output$formula <- renderPrint({
    fml_maker(input$class, input$C, input$db)
  })
  
  output$mzvals <- renderPrint({
    fml <- fml_maker(input$class, input$C, input$db)
    mz_calculator(input$class, fml)
  })
  
  output$cmps_pred <- renderPrint({
    cmps_db$name[unlist(matchWithPpm(input$mzs, cmps_db[,input$pol], ppm = input$ppms))]
  })
  
  # MSMS POS ----
  # to add mPA, dmPA + from PC to DG + from TG;O
  output$ms2_thr_pos <- renderPlot({
    fml <- fml_maker(input$class, input$C, input$db)
    mz <- mz_calculator(input$class, fml)
    mz <- mz[grep("]\\+", colnames(mz))]
    if(input$class == "CAR"){
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements(fml, "C3H9N")), "[M+H]+")),
        i = 100,
        ad = "[M+H-C3H9N]+"
      )
    } else if(input$class %in% c("HexCer", "HexCer;O3", "HexCer;O4")){ ## HexCer ----
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O5")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H12O6")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H12O6H2O")), "[M+H]+")),
               as.numeric(calculateMass(subtractElements(addElements(sn$formula[idx1], "NH2"), "O"))),
               as.numeric(calculateMass(subtractElements(addElements(sn$formula[idx2], "N"), "O2"))),
               #as.numeric(calculateMass(subtractElements(addElements(sn$formula[idx2], "NH2"), "O"))),
               as.numeric(calculateMass(addElements(sn$formula[idx2], "C3H4N"))),
               as.numeric(calculateMass(addElements(sn$formula[idx2], "C3H6NO")))
        ),        
        i = c(25, 70, 100, 12, 30, #2, 
              5, 15),
        ad = c("[M+H-hexose]+", "[M+H-hexose-H2O]+", "[M+H-hexose-2(H2O)]+",
               paste0("[", sn$sn[idx1], "+NH2-O]+"),
               paste0("[", sn$sn[idx2], "+N-O2]+"),
               #paste0("[", sn$sn[idx2], "+NH2-O]+"),
               paste0("[", sn$sn[idx2], "+C3H4N]+"),
               paste0("[", sn$sn[idx2], "+C3H6NO]+")
        )
      )
    } else if(input$class == "LPA"){ ## LPA ----
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))),
        i = c(100, 20),
        ad = c("[M+H-H2O]+", "[M+H]+")
      )
    }else if(input$class == "LPC"){
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements(fml, "H2O")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass("C5H14NO4P")), "[M+H]+")),
        i = c(100, 50),
        ad = c("[M+H-H2O]+", "[PC+H]+")
      )
    } else if(input$class == "LPE"){
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements(fml, "H2O")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C2H8NO4P"))), "[M+H]+")),
        i = c(100, 50),
        ad = c("[M+H-H2O]+", "[M+H-PE]+")
      )
    } else if(input$class == "LPS"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C3H7NO3"))), "[M+H]+"),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C3H7NO3HPO3"))), "[M+H]+")
        ),
        i = c(100, 20, 30),
        ad = c("[M+H-H2O]+", "[M+H-S]+", "[M+H-PS]+")
      )
    }else if(input$class == "PA"){
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "H3PO4"))), "[M+H]+")),
        i = c(30, 100),
        ad = c("[M+H]+", "[M+H-PA]+")
      )
      # if input$C >= 40....
    } else if(input$class == "PC"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mz - sn$mass[idx1]),
               as.numeric(mz - sn$mass[idx2]),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O"))
        ), 
        i = c(90, 40, 10, 100),
        ad = c(paste0("[M+H-", input$sn1, "-H2O]+"),
               paste0("[M+H-", input$sn2, "-H2O]+"),
               paste0("[M+H-", input$sn1, "]+"),
               paste0("[M+H-", input$sn2, "]+")
        ))
      if(input$sn1 == input$sn2){
        sps <- sps[c(2, 4),]
      }
    } else if(input$class == "PE"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements(fml, "C2H8NO4P"))), "[M+H]+"),
        i = 100,
        ad = "[M+H-PE]+"
      )
    } else if(input$class == "PG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      sps <- data.frame(
        mz = c(mz,
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C3H8O3HPO3")), "[M+H]+"))
        ),
        i = c(60, 100),
        ad = c("[M+H]+", "[M+H-PG]+")
      )
    } else if(input$class == "PS"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements(fml, "C3H7NO3HPO3"))), "[M+H]+"),
        i = 100,
        ad = "[M+H-PS]+"
      )
    } else if(input$class == "MGDG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(mz,
               as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O5")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H12O6")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx1], "C3H4O")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx2], "C3H4O")), "[M+H]+"))
        ),
        i = c(10, 20, 100, 30, 15, 5),
        ad = c("[M+H]+", "[M+H-H2O]+", "[M+H-hexose]+", "[M+H-hexose-H2O]+",
               paste0("[", input$sn1, "+H+C3H4O]+"),
               paste0("[", input$sn2, "+H+C3H4O]+"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[-nrow(sps),]
      }
    } else if(input$class == "DGDG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(mz,
               as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O5C6H10O5")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O5C6H12O6")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx1], "C3H4O")), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx2], "C3H4O")), "[M+H]+"))
        ),
        i = c(10, 20, 100, 30, 15, 5),
        ad = c("[M+H]+", "[M+H-H2O]+", "[M+H-2(hexose)]+", "[M+H-2(hexose)-H2O]+",
               paste0("[", input$sn1, "+H+C3H4O]+"),
               paste0("[", input$sn2, "+H+C3H4O]+"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[-nrow(sps),]
      }
    } else if(input$class == "DGGA"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+NH4]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(#mz,
          as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
          as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O7")), "[M+H]+")),
          as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H8O6")), "[M+H]+")),
          as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx1], "C3H4O")), "[M+H]+")),
          as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx2], "C3H4O")), "[M+H]+"))
        ),
        i = c(10, 100, 20, 80, 40),
        ad = c("[M+H-H2O]+", "[M+H-glucuronide-H2O]+", "[M+H-glucuronide]+",
               paste0("[", input$sn1, "+H+C3H4O]+"),
               paste0("[", input$sn2, "+H+C3H4O]+")
        )
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(1, 2, 4),]
      }
    } else if(input$class == "DGTS"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+H]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mz - sn$mass[idx1]),
               as.numeric(mz - sn$mass[idx2]),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O")),
               as.numeric(mass2mz(calculateMass(subtractElements(addElements(addElements("C4H9NO3", "CH3CH3CH3"), "C3H8O3"), "H2OH3"))))
        ), 
        i = c(60, 100, 40, 90, 25),
        ad = c(paste0("[M+H-", input$sn1, "-H2O]+"),
               paste0("[M+H-", input$sn2, "-H2O]+"),
               paste0("[M+H-", input$sn1, "]+"),
               paste0("[M+H-", input$sn2, "]+"),
               "[GTS+H]+"
        ))
      if(input$sn1 == input$sn2){
        sps <- sps[c(2, 4),]
      }
    } else if(input$class == "MG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+NH4]+"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(fml), "[M+H]+")), 
        i = 100,
        ad = "[M+H]+"
      )
    } else if(input$class == "DG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+NH4]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M+H]+")), 
               as.numeric(mass2mz(calculateMass(fml), "[M+H-H2O]+")),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+") - sn$mass[idx1]),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+") - sn$mass[idx2])
        ), 
        i = c(30, 100, 60, 40),
        ad = c("[M+H]+",
               "[M+H-H2O]+",
               paste0("[M+H-", input$sn1, "-H2O]+"),
               paste0("[M+H-", input$sn2, "-H2O]+")
        ))
      if(input$sn1 == input$sn2){
        sps <- sps[-nrow(sps),]
      } 
    } else if(input$class == "TG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+NH4]+"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      idx3 <- which(sn$sn == input$sn3)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M+H]+")),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+") - sn$mass[idx1]),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+") - sn$mass[idx2]),
               as.numeric(mass2mz(calculateMass(fml), "[M+H]+") - sn$mass[idx3])
        ), 
        i = c(30, 100, 80, 60),
        ad = c("[M+H]+",
               paste0("[M+H-", input$sn1, "-H2O]+"),
               paste0("[M+H-", input$sn2, "-H2O]+"),
               paste0("[M+H-", input$sn3, "-H2O]+")
        ))
      if(input$sn1 == input$sn2 & input$sn2 == input$sn3){
        sps <- sps[c(1, 2),]
      } else if(input$sn2 == input$sn3){
        sps <- sps[c(1, 2, 3),]
      }else if(input$sn1 == input$sn2){
        sps <- sps[c(1, 2, 4),]
      }
    }
    if(exists("sps")){
      sps$mz <- as.numeric(sps$mz)
      sps$i <- as.numeric(sps$i)
      plot(sps$mz, sps$i, type = "h", bty = "l", 
           xlim = c(50, mz + 2), ylim = c(0, 110), xlab = "m/z", ylab = "intensity", 
           main = paste("Precursor m/z", sprintf("%.4f", round(mz, 4))))
      text(sps$mz, sps$i, paste(sprintf("%.4f", sps$mz), "\n", sps$ad), pos = 3)
    }
  })
  
  # MSMS NEG ----
  
  output$ms2_thr_neg <- renderPlot({
    fml <- fml_maker(input$class, input$C, input$db)
    mz <- mz_calculator(input$class, fml)
    mz <- mz[grep("]\\-", colnames(mz))]
    if(input$class == "pHexFA"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      idx1 <- which(sn$C == input$C & sn$db == input$db)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx1], "[M-H]-"))),
        i = c(100, 30),
        ad = c("[M-H]-", paste0("[", sn$sn[idx1],"-H]-"))
      )
    } else if(input$class == "SM"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      idx1 <- which(sn$sn == input$sn1)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH3")), "[M]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C5H13NO")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, subtractElements(addElements(sn$formula[idx1], "N"), "H2O"))), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(subtractElements(fml, subtractElements(sn$formula[idx1], "H2O")), "C5H13NO")), "[M-H]-"))
        ),        
        i = c(100, 20, 100, 10),
        ad = c("[M-CH3]-", "[M-H-C5H13NO]-", paste0("[M-H-", sn$sn[idx1],"]-"),
               paste0("[M-H-C5H13NO-", sn$sn[idx1],"]-"))
      )
    } else if(input$class %in% c("Cer", "Cer;O3", "Cer;O4")){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "H2O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH2O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH3OH")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH2OH2O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx1], "C2NH3")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(addElements(sn$formula[idx1], "NH3C2"), "O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(addElements(sn$formula[idx1], "NH3C2"), "CO")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(sn$formula[idx2], "H2O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(sn$formula[idx1]), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(addElements(subtractElements(sn$formula[idx1], "O"), "NH")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(sn$formula[idx2], "C2H4O")), "[M-H]-"))
        ),
        i = c(80, 5, 100, 40, 15, 20, 75, 5, 60, 40, 25, 30),
        ad = c("[M-H]-", "[M-H-H2O]-", "[M-H-CH2O]-", "[M-H-CH3OH]-", "[M-H-CH2O-H2O]-", 
               paste("S", sn$sn[idx1]), paste("T", sn$sn[idx1]), paste("X", sn$sn[idx1]), 
               paste("R", sn$sn[idx2]), paste0("[", sn$sn[idx1], "-H]-"), paste("U", sn$sn[idx1]), paste("P", sn$sn[idx2]))
      )
    } else if(input$class %in% c("HexCer", "HexCer;O3", "HexCer;O4")){ ## HexCer ----
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H10O5")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H12O6")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(addElements(sn$formula[idx1], "C2H3N")), "[M-H]-"))
        ),        
        i = c(100, 50, 100, 20),
        ad = c("[M-H]-", "[M-H-hexose]-", "[M-H-hexose-H2O]-",
               paste0("[", sn$sn[idx1], "-H+C2H3N]-"))
      )
    } else if(input$class == "LPA"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements("C3H9O6P", "H2O")), "[M-H]-")),        
        i = 100,
        ad = "[GP-H-H2O]-"
      )
    } else if(input$class == "LPC" | input$class == "PC"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH3")), "[M]-")),        
        i = 100,
        ad = "[M-CH3]-"
      )
    } else if(input$class == "LPE"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(
          paste0("C", input$C, "H", input$C*2 - (2*input$db), "O2")), "[M-H]-")), 
        i = 100,
        ad = paste0("[", input$C, ":", input$db, "-H]-")
      )
    } else if(input$class == "LPG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$C == input$C & sn$db == input$db)
      sps <- data.frame(
        mz = as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")), 
        i = 100,
        ad = paste0("[", sn$sn[idx1], "-H]-")
      )
    } else if(input$class == "LPI"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$C == input$C & sn$db == input$db)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements("C6H13O9P", "H2O")), "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mz - sn$mass[idx1]),
               as.numeric(mass2mz(calculateMass(subtractElements(fml, "C6H12O6")), "[M-H]-"))
        ),        
        i = c(40, 100, 40, 50),
        ad = c("[IP-H-H2O]-",
               paste0("[", sn$sn[idx1], "-H]-"),
               paste0("[M-H-", sn$sn[idx1], "-H2O]-"),
               "[M-H-hexose-H2O]-")
      )
    } else if(input$class == "LPS"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$C == input$C & sn$db == input$db)
      sps <- data.frame(
        mz = as.numeric(mz - calculateMass(subtractElements("C3H7NO3", "H2O"))), 
        i = 100,
        ad = "[M-H-Ser]-"
      )
    } else if(input$class == "PA"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-")),
               as.numeric(mz - sn$mass[idx1]),
               as.numeric(mz - sn$mass[idx2]),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O"))
        ), 
        i = c(80, 40, 40, 100, 5, 60),
        ad = c(paste0("[", input$sn1, "-H]-"), 
               paste0("[", input$sn2, "-H]-"),
               paste0("[M-H-", input$sn1, "-H2O]-"),
               paste0("[M-H-", input$sn2, "-H2O]-"),
               paste0("[M-H-", input$sn1, "]-"),
               paste0("[M-H-", input$sn2, "]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(2,4,6),]
      }
    } else if(input$class == "PE"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-")),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O"))
        ), 
        i = c(50, 100, 10, 20),
        ad = c(paste0("[", input$sn1, "-H]-"), 
               paste0("[", input$sn2, "-H]-"),
               paste0("[M-H-", input$sn1, "]-"),
               paste0("[M-H-", input$sn2, "]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(2,4),]
      }
    } else if(input$class == "PG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-")),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O") - calculateMass("C3H8O3")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O") - calculateMass("C3H8O3"))
        ), 
        i = c(50, 100, 5, 20, 5, 20),
        ad = c(paste0("[", input$sn1, "-H]-"), 
               paste0("[", input$sn2, "-H]-"),
               paste0("[M-H-", input$sn1, "]-"),
               paste0("[M-H-", input$sn2, "]-"),
               paste0("[M-H-", input$sn1, "-glycerol]-"),
               paste0("[M-H-", input$sn2, "-glycerol]-")
        )
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(2,4,6),]
      }
    } else if(input$class == "PI"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-")),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O")),
               as.numeric(mz - sn$mass[idx1] + calculateMass("H2O") - calculateMass("C6H12O6")),
               as.numeric(mz - sn$mass[idx2] + calculateMass("H2O") - calculateMass("C6H12O6")),
               as.numeric(mz - sn$mass[idx1]),
               as.numeric(mz - sn$mass[idx2])
        ), 
        i = c(50, 40, 5, 20, 10, 50, 15, 100),
        ad = c(paste0("[", input$sn1, "-H]-"), 
               paste0("[", input$sn2, "-H]-"),
               paste0("[M-H-", input$sn1, "]-"),
               paste0("[M-H-", input$sn2, "]-"),
               paste0("[M-H-", input$sn1, "-inositol]-"),
               paste0("[M-H-", input$sn2, "-inositol]-"),
               paste0("[M-H-", input$sn1, "-H2O]-"),
               paste0("[M-H-", input$sn2, "-H2O]-")
        )
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(2,4,6,8),]
      }
    } else if(input$class == "PS"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mz - calculateMass(subtractElements("C3H7NO3", "H2O"))),
               as.numeric(mass2mz(calculateMass(subtractElements(subtractElements(fml, "C3H5NO2"), sn$formula[idx1])), "[M-H]-")), 
               as.numeric(mass2mz(calculateMass(subtractElements(subtractElements(fml, "C3H5NO2"), sn$formula[idx2])), "[M-H]-"))), 
        i = c(100, 20, 10),
        ad = c("[M-H-Ser]-", paste0("[M-H-Ser-", sn$sn[idx1], "]-"), paste0("[M-H-Ser-", sn$sn[idx2], "]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(1:2),]
      }
    } else if(input$class == "MGDG"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M+CHO2]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(fml), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(subtractElements(sn$formula[idx1], "H2O"))),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(subtractElements(sn$formula[idx2], "H2O"))),
               as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-"))
        ),
        i = c(100, 20, 10, 30, 60),
        ad = c("[M-H]-", paste0("[M-H-", sn$sn[idx1],"]-"), paste0("[M-H-", sn$sn[idx2],"]-"), 
               paste0("[", sn$sn[idx1],"-H]-"), paste0("[", sn$sn[idx2], "-H]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(1,2,4),]
      }
    } else if(input$class == "DGGA"){
      fml <- fml_maker(input$class, input$C, input$db)
      mz <- as.numeric(mass2mz(calculateMass(fml), "[M-H]-"))
      idx1 <- which(sn$sn == input$sn1)
      idx2 <- which(sn$sn == input$sn2)
      sps <- data.frame(
        mz = c(as.numeric(mass2mz(calculateMass(subtractElements(fml, "H2O")), "[M-H]-")),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(subtractElements(sn$formula[idx1], "H2O"))),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(subtractElements(sn$formula[idx2], "H2O"))),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(sn$formula[idx1])),
               as.numeric(mass2mz(calculateMass(fml), "[M-H]-") - calculateMass(sn$formula[idx2])),
               as.numeric(mass2mz(sn$mass[idx1], "[M-H]-")),
               as.numeric(mass2mz(sn$mass[idx2], "[M-H]-"))
        ),
        i = c(20, 35, 100, 5, 20, 80, 40),
        ad = c("[M-H-H2O]-", 
               paste0("[M-H-", sn$sn[idx1],"]-"), paste0("[M-H-", sn$sn[idx2],"]-"), 
               paste0("[M-H-", sn$sn[idx1],"-H2O]-"), paste0("[M-H-", sn$sn[idx2],"-H2O]-"), 
               paste0("[", sn$sn[idx1],"-H]-"), paste0("[", sn$sn[idx2], "-H]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(1,3,5,6),]
      }
    } 
    if(exists("sps")){
      sps$mz <- as.numeric(sps$mz)
      sps$i <- as.numeric(sps$i)
      plot(sps$mz, sps$i, type = "h", bty = "l", 
           xlim = c(50, mz + 2), ylim = c(0, 110), xlab = "m/z", ylab = "intensity", 
           main = paste("Precursor m/z", sprintf("%.4f", round(mz, 4))))
      text(sps$mz, sps$i, paste(sprintf("%.4f", sps$mz), "\n", sps$ad), pos = 3)
    }
  })
  
  output$kmd <- renderPlotly({
    cmps <- cmps[cmps$class %in% input$kmd_class,]
    
    idx <- which(cmps$standard == TRUE)
    plot_ly(data = cmps) %>%
      add_markers(data = cmps[idx, ], x = ~rtime, y = ~KMD, size = 3,
                  color = cmps[idx, input$kmd_color],
                  symbol = as.character(cmps$standard[idx] == TRUE), text = NULL,
                  symbols = 4, mode = "markers", showlegend = FALSE, 
                  text = paste("mz: ", round(cmps$mass[idx], 4), 
                               "\n KMD: ", round(cmps$KMD[idx], 8),
                               "\n RT: ", round(cmps$rtime[idx], 2),
                               "\n Abbreviation:", 
                               paste0(cmps$class[idx], cmps$C[idx], ":", cmps$db[idx]),
                               "\n name: ", cmps$compound[idx]
                  ), 
                  hoverinfo = "text") %>%
      add_markers(data = cmps, x = ~rtime, y = ~KMD, 
                  text = paste("mz: ", round(cmps$mass, 4), 
                               "\n KMD: ", round(cmps$KMD, 8),
                               "\n RT: ", round(cmps$rtime, 2),
                               "\n Abbreviation:", 
                               "\n", paste0(cmps$class, cmps$C, ":", cmps$db),
                               "\n name: ", cmps$compound
                  ), 
                  hoverinfo = "text", size = 2,
                  color = cmps[,input$kmd_color])
  })
  
  # ESI+ ----
  output$posfrag1add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag1) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag1, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  output$posfrag2add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag2) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag2, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  output$posfrag3add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag3) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag3, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  output$posfrag4add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag4) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag4, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  output$posfrag5add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag5) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag5, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  output$posfrag6add <- renderPrint({
    idx1 <- c(which(abs((input$posprec - input$posfrag6) - mzdif.pos1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$posfrag6, mzdif.pos2$dif, ppm = 10)))
    unlist(c(mzdif.pos1$add[idx1], mzdif.pos2[idx2]))
  })
  
  
  # ESI- -----
  output$negfrag1add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag1) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag1, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  output$negfrag2add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag2) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag2, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  output$negfrag3add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag3) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag3, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  output$negfrag4add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag4) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag4, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  output$negfrag5add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag5) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag5, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  output$negfrag6add <- renderPrint({
    idx1 <- c(which(abs((input$negprec - input$negfrag6) - mzdif.neg1$dif) < 0.01))
    idx2 <- c(unlist(matchWithPpm(input$negfrag6, mzdif.neg2$dif, ppm = 10)))
    unlist(c(mzdif.neg1$add[idx1], mzdif.neg2[idx2]))
  })
  
  
  
  
}

shinyApp(ui = ui, server = server)
options(repos = BiocManager::repositories())
options("repos")

cmps <- readxl::read_xlsx("compounds.xlsx")

library(shiny)
library(tidyverse)
library(MetaboCoreUtils)
library(plotly)


.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

# to add mPA, dmPA
fml_maker <- function(class, C, db){
  case_when(
    # 1 FA chain -----------------------------------------------------------
    class == "FA" ~ paste0("C", C, "H", C*2 - (2*db), "O2"),
    class == "CAR" ~ paste0("C", C + 7, "H", C*2 - (2*db) + 13, "NO4"),
    class == "LPA" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 7, "O7P"),
    class == "LPC" ~ paste0("C", C + 8, "H", C*2 - (2*db) + 18, "NO7P"),
    class == "LPE" ~ paste0("C", C + 5, "H", C*2 - (2*db) + 12, "NO7P"),
    class == "LPG" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 13, "O9P"),
    class == "LPI" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 17, "O12P"),
    class == "LPS" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 12, "NO9P"),
    class == "MG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  6, "O4"),
    
    # 2 FA chain -----------------------------------------------------------
    class == "PA"  ~ paste0("C", C + 3, "H", C*2 - (2*db) +  5, "O8P"),
    class == "PC"  ~ paste0("C", C + 8, "H", C*2 - (2*db) + 16, "NO8P"),
    class == "PE"  ~ paste0("C", C + 5, "H", C*2 - (2*db) + 10, "NO8P"),
    class == "PG"  ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "O10P"),
    class == "PI"  ~ paste0("C", C + 9, "H", C*2 - (2*db) + 15, "O13P"),
    class == "PS"  ~ paste0("C", C + 6, "H", C*2 - (2*db) + 10, "NO10P"),
    
    class == "SM"       ~ paste0("C", C + 4, "H", C*2 - (2*db) + 11,"N2O6P"),
    class == "Cer"      ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO3"),
    class == "Cer;O3"   ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO4"),
    class == "Cer;O4"   ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO5"),
    class == "HexCer"   ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO8"),
    class == "LactCer"  ~ paste0("C", C +12, "H", C*2 - (2*db) + 21, "NO13"),
    
    class == "MGDG" ~ paste0("C", C +  9, "H", C*2 - (2*db) + 14, "O10"),
    class == "DGDG" ~ paste0("C", C + 15, "H", C*2 - (2*db) + 24, "O15"),
    class == "DGTS" ~ paste0("C", C + 10, "H", C*2 - (2*db) + 17, "NO7"),
    
    class == "DG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  4, "O5"),
    
    # 3 FA chain -----------------------------------------------------------
    class == "TG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O6")
  )
}



cmps_db <- c()

cls <- c("FA", "CAR", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "MG")
C <- seq(from = 12, to = 24, by = 1)
db <- seq(from = 0, to = 6, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}

cls <- c("SM", "Cer", "Cer;O3", "Cer;O4", "HexCer", "LactCer", 
         "PA", "PC", "PE", "PG", "PI", "PS", "MGDG", "DGDG", "DGTS", "DG")
C <- seq(from = 12*2, to = 28*2, by = 1)
db <- seq(from = 0, to = 6*2, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}

cls <- "TG"
C <- seq(from = 12*3, to = 24*3, by = 1)
db <- seq(from = 0, to = 6*3, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}


cmps_db <- data.frame(
  name = cmps_db,
  class = gsub("\\ .*", "", cmps_db),
  C = gsub(":.*", "", gsub(".*\\ ", "", cmps_db)),
  db = gsub(".*:", "", gsub(".*\\ ", "", cmps_db))
)
cmps_db$C <- as.numeric(cmps_db$C)
cmps_db$db <- as.numeric(cmps_db$db)
cmps_db$formula <- fml_maker(cmps_db$class, cmps_db$C, cmps_db$db)
cmps_db$mass <- calculateMass(cmps_db$formula)
cmps_db$pos <- cmps_db$neg <- NA
idx <- which(cmps_db$class %in% c("FA", "CAR", 
                                  "SM", "Cer", "Cer;O3", "Cer;O4", "HexCer", "LactCer", 
                                  "LPC", "LPE", "LPS", "PC", "PE", "PS",
                                  "DGTS"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+H]+")
idx <- which(cmps_db$class %in% c("LPA", "LPG", "LPI",
                                  "PA", "PG", "PI", "MGDG", "DGDG", 
                                  "MG", "DG", "TG"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+NH4]+")
idx <- which(cmps_db$class %in% c("FA", "LPA", "LPE", "LPG", "LPI", "LPS", 
                                  "PA", "PE", "PG", "PI", "PS", 
                                  "MG", "DG", "TG"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M-H]-")
idx <- which(cmps_db$class %in% c("CAR", 
                                  "SM", "Cer", "Cer;O3", "Cer;O4", "HexCer", "LactCer", 
                                  "LPC", "PC", "MGDG", "DGDG", "DGTS"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M+CHO2]-")



mz_calculator <- function(class, fml){
  if(class == "FA"){
    mass2mz(calculateMass(fml), c("[M-H]-"))
  } else if(class %in% c("Cer", "Cer;O3", "Cer;O4", "HexCer", "LactCer")){
    mass2mz(calculateMass(fml), c("[M+CHO2]-"))
  } else if(class %in% c("MG", "DG", "TG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+"))
  } else if(class %in% c("LPE", "LPS", "PE", "PS")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M-H]-"))
  } else if(class %in% c("SM", "LPC", "PC", "DGTS")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M+CHO2]-"))
  } else if(class %in% c("CAR", "LPA", "LPG", "LPI", "PA", "PG", "PI")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M-H]-"))
  } else if(class %in% c("MGDG", "DGDG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M+CHO2]-"))
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
  "C" = rep(seq(12, 24), each = 4),
  "db" = rep(seq(0, 3), 13)
)
sn$formula <- paste0("C", sn$C, "H", sn$C*2 - 2*sn$db, "O2")
sn$sn <- paste0(sn$C, ":", sn$db)
sn$mass <- NA
sn.list <- list()
for(i in seq(nrow(sn))){
  sn$mass[i] <- calculateMass(sn$formula[i])
  sn.list[[i]] <- sn$sn[i]
  names(sn.list)[[i]] <- sn$sn[i]
}

mzdif.pos <- data.frame(rbind(
  c(calculateMass("NH3"), "loss NH3 -> LPA / PA / PG / MGDG / MG / DG / TG"),
  c(calculateMass("NH3H2O"), "loss NH3 & H2O -> LPA / MGDG / DG"),
  c(calculateMass("H3PO4NH3"), "loss NH3 & PA -> PA"),
  c(calculateMass("CH3"), "loss CH3"),
  cbind(sn$mass + calculateMass("NH3"), paste("loss NH3 &", sn$sn, "-> DG / TG")),
  c(calculateMass("H2O"), "loss H2O -> LPC / LPS"),
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
  
  cbind(sn$mass, paste0("loss '", sn$sn, "-H2O' -> PC / DGTS")),
  cbind(sn$mass - calculateMass("H2O"), paste("loss ", sn$sn, " -> PC/ DGTS")),
  cbind(mass2mz(calculateMass(addElements(sn$formula, "C3H4O")), "[M+H]+"), paste0("[", sn$sn, "+H+C3H4O]+ -> MGDG / DGDG"))
))
colnames(mzdif.pos) <- c("dif", "add")
mzdif.pos$dif <- as.numeric(mzdif.pos$dif)

mzdif.neg <- data.frame(rbind(
  c(calculateMass("HCOOH"), "loss HCOOH -> CER / MGDG / DGDG")
))
colnames(mzdif.neg) <- c("dif", "add")
mzdif.neg$dif <- as.numeric(mzdif.neg$dif)

# UI ------------------------------------------------------------------------
ui <- navbarPage(
  "Lipidomics",
  
  tabPanel(
    "Main Panel",
    column(6, h3("Formula"),
           fluidRow(
             column(2, selectInput(inputId = "class", label = "Lipid class:",
                                   choices = list("FA" = "FA",
                                                  "CAR" = "CAR", 
                                                  "SM" = "SM",
                                                  "Cer" = "Cer",
                                                  "Cer;O3" = "Cer;O3",
                                                  "Cer;O4" = "Cer;O4",
                                                  "HexCer" = "HexCer",
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
  tabPanel(
    "Kendrick Mass Defect (KMD)",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("kmd_class", label = "Lipid classes:", 
                           choices = list("FA" = "FA",
                                          "CAR" = "CAR", 
                                          "SM" = "SM",
                                          "Cer" = "Cer",
                                          "Cer;O3" = "Cer;O3",
                                          "Cer;O4" = "Cer;O4",
                                          "HexCer" = "HexCer",
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
                                          "DGTS" = "DGTS",
                                          "MG" = "MG",
                                          "DG" = "DG",
                                          "TG" = "TG"),
                           selected = "TG"),
        radioButtons("kmd_color", label = "Color by:",
                     choices = list("Class" = "class", 
                                    "Double bonds" = "db"),
                     selected = "class")
      ),
      mainPanel(
        plotlyOutput("kmd")
      )
    )), # close tabPanel "KMD"
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
    } else if(input$class == "LPA"){
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
               as.numeric(mz - sn$mass[idx1]),
               as.numeric(mz - sn$mass[idx2])
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
  
  output$ms2_thr_neg <- renderPlot({
    fml <- fml_maker(input$class, input$C, input$db)
    mz <- mz_calculator(input$class, fml)
    mz <- mz[grep("]\\-", colnames(mz))]
    if(input$class == "LPC"){
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(subtractElements(fml, "CH3")), "[M]-")),        
        i = 100,
        ad = "[M-CH3]-"
      )
    } else if(input$class == "LPE"){
      sps <- data.frame(
        mz = as.numeric(mass2mz(calculateMass(
          paste0("C", input$C, "H", input$C*2 - (2*input$db), "O2")), "[M-H]-")), 
        i = 100,
        ad = paste0("[", input$C, ":", input$db, "H]-")
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
        ad = c(paste0("[", input$sn1, "H]-"), 
               paste0("[", input$sn2, "H]-"),
               paste0("[M-H-", input$sn1, "-H2O]-"),
               paste0("[M-H-", input$sn2, "-H2O]-"),
               paste0("[M-H-", input$sn1, "]-"),
               paste0("[M-H-", input$sn2, "]-"))
      )
      if(input$sn1 == input$sn2){
        sps <- sps[c(2,4,6),]
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
  
  output$posfrag5add <- renderPrint({
    idx <- c(which(abs((input$posprec - input$posfrag5) - mzdif.pos$dif) < 0.01),
             unlist(matchWithPpm(input$posfrag5, mzdif.pos$dif, ppm = 10)))
    mzdif.pos$add[idx]
  })
  
  output$posfrag6add <- renderPrint({
    idx <- c(which(abs((input$posprec - input$posfrag6) - mzdif.pos$dif) < 0.01),
             unlist(matchWithPpm(input$posfrag6, mzdif.pos$dif, ppm = 10)))
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
  
}

shinyApp(ui = ui, server = server)
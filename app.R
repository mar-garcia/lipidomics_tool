options(repos = BiocManager::repositories())
#options("repos")


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


fml_maker <- function(class, C, db){
  case_when(
    # 1 FA chain -----------------------------------------------------------
    class == "FA" ~ paste0("C", C, "H", C*2 - (2*db), "O2"),
    class == "CAR" ~ paste0("C", C + 7, "H", C*2 - (2*db) + 13, "NO4"),
    class == "LPA" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 7, "O7P"),
    class == "LPC" ~ paste0("C", C + 8, "H", C*2 - (2*db) + 18, "NO7P"),
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
    class == "HexCer"   ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO8"),
    class == "LactCer"  ~ paste0("C", C +12, "H", C*2 - (2*db) + 21, "NO13"),
    
    class == "MGDG" ~ paste0("C", C +  9, "H", C*2 - (2*db) + 14, "O10"),
    class == "DGDG" ~ paste0("C", C + 15, "H", C*2 - (2*db) + 24, "O15"),
    
    class == "DG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  4, "O5"),
    
    # 3 FA chain -----------------------------------------------------------
    class == "TG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O6")
  )
}



cmps_db <- c()

cls <- c("FA", "CAR", "LPA", "LPC", "LPG", "LPI", "LPS", "MG")
C <- seq(from = 12, to = 24, by = 1)
db <- seq(from = 0, to = 6, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}

cls <- c("SM", "Cer", "HexCer", "LactCer", 
         "PA", "PC", "PE", "PG", "PI", "PS", "MGDG", "DGDG", "DG")
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
idx <- which(cmps_db$class %in% c("FA", "CAR", "SM", "Cer", "HexCer", "LactCer", 
                                  "LPC", "LPS", "PC", "PE", "PS"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+H]+")
idx <- which(cmps$class %in% c("LPA", "LPG", "LPI",
                               "PA", "PG", "PI", "MGDG", "DGDG", 
                               "MG", "DG", "TG"))
cmps_db$pos[idx] <- mass2mz(cmps_db$mass[idx], "[M+NH4]+")
idx <- which(cmps_db$class %in% c("FA", "LPA", "LPG", "LPI", "LPS", 
                                  "PA", "PE", "PG", "PI", "PS", 
                                  "MG", "DG", "TG"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M-H]-")
idx <- which(cmps_db$class %in% c("SM", "CAR", "Cer", "HexCer", "LactCer", 
                                  "LPC", "PC", "MGDG", "DGDG"))
cmps_db$neg[idx] <- mass2mz(cmps_db$mass[idx], "[M+CHO2]-")



mz_calculator <- function(class, fml){
  if(class == "FA"){
    mass2mz(calculateMass(fml), c("[M-H]-"))
  } else if(class %in% c("Cer", "HexCer", "LactCer")){
    mass2mz(calculateMass(fml), c("[M+CHO2]-"))
  } else if(class %in% c("MG", "DG", "TG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+"))
  } else if(class %in% c("LPS", "PE", "PS")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M-H]-"))
  } else if(class %in% c("SM", "LPC", "PC")){
    mass2mz(calculateMass(fml), c("[M+H]+", "[M+CHO2]-"))
  } else if(class %in% c("CAR", "LPA", "LPG", "LPI", "PA", "PG", "PI")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M-H]-"))
  } else if(class %in% c("MGDG", "DGDG")){
    mass2mz(calculateMass(fml), c("[M+NH4]+", "[M+CHO2]-"))
  }
}

ref <- c(12, 1, 16, 14, 31, 32)
names(ref) <- c("C", "H", "O", "N", "P", "S")


cmps <- read.csv("compounds.csv")
cmps <- cmps[cmps$class %in% c("FA", "CAR", "SM", "Cer", "HexCer", "LactCer", 
                               "LPA", "LPC", "LPG", "LPI", "LPS", 
                               "PA", "PC", "PE", "PG", "PI", "PS",
                               "MGDG", "DGDG", "MG", "DG", "TG"),]

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


# UI ------------------------------------------------------------------------
ui <- navbarPage(
  "Lipidomics",
  
  tabPanel(
    "Panel 1",
    column(6, h3("Formula"),
           fluidRow(
             column(2, selectInput(inputId = "class", label = "Lipid class:",
                                   choices = list("FA" = "FA",
                                                  "CAR" = "CAR", 
                                                  "SM" = "SM",
                                                  "Cer" = "Cer",
                                                  "HexCer" = "HexCer",
                                                  "LactCer" = "LactCer",
                                                  "LPA" = "LPA",
                                                  "LPC" = "LPC",
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
           ) # close right column
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
                                          "HexCer" = "HexCer",
                                          "LactCer" = "LactCer",
                                          "LPA" = "LPA",
                                          "LPC" = "LPC",
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
                                          "MG" = "MG",
                                          "DG" = "DG",
                                          "TG" = "TG"),
                           selected = c("FA", "CAR", 
                                        "SM", "Cer", "HexCer", "LactCer", 
                                        "LPA", "LPC", "LPG", "LPI", "LPS",
                                        "PA", "PC", "PE", "PG", "PI", "PS",
                                        "MGDG", "DGDG", "MG", "DG", "TG")),
        radioButtons("kmd_color", label = "Color by:",
                     choices = list("Class" = "class", 
                                    "Double bonds" = "db"),
                     selected = "class")
      ),
      mainPanel(
        plotlyOutput("kmd")
      )
    )) # close tabPanel "KMD"
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
  
}

shinyApp(ui = ui, server = server)
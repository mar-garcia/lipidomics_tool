fml_maker <- function(class, C = 0, db = 0, M = 0){
  case_when(
    # 1 FA chain -----------------------------------------------------------
    class == "FA" ~ paste0("C", C, "H", C*2 - (2*db), "O2"),
    class == "FA;O" ~ paste0("C", C, "H", C*2 - (2*db), "O3"),
    class == "FA;O2" ~ paste0("C", C, "H", C*2 - (2*db), "O4"),
    class == "FA;O3" ~ paste0("C", C, "H", C*2 - (2*db), "O5"),
    class == "FA;O4" ~ paste0("C", C, "H", C*2 - (2*db), "O6"),
    class == "pHexFA" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 12, "O7"),
    class == "pPentFA"~ paste0("C", C + 5, "H", C*2 - (2*db) + 10, "O6"),
    class == "CAR" ~ paste0("C", C + 7, "H", C*2 - (2*db) + 13, "NO4"),
    class == "LPA" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 7, "O7P"),
    class == "LPC" ~ paste0("C", C + 8, "H", C*2 - (2*db) + 18, "NO7P"),
    class == "LPE" ~ paste0("C", C + 5, "H", C*2 - (2*db) + 12, "NO7P"),
    class == "LPG" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 13, "O9P"),
    class == "LPI" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 17, "O12P"),
    class == "LPS" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 12, "NO9P"),
    class == "DGMG"~ paste0("C", C +15, "H", C*2 - (2*db) + 26, "O14"),
    class == "MG"  ~ paste0("C", C + 3, "H", C*2 - (2*db) +  6, "O4"),
    class == "ST" ~ paste0("C", C +  29, "H", C*2 - (2*db) + 48, "O2"),
    class == "Glc-ST" ~ paste0("C", C +  35, "H", C*2 - (2*db) + 58, "O7"),
    
    # 2 FA chain -----------------------------------------------------------
    class == "PA"  ~ paste0("C", C + 3, "H", C*2 - (2*db) +  5, "O8P"),
    class == "mPA" ~ paste0("C", C + 4, "H", C*2 - (2*db) +  7, "O8P"),
    class == "dmPA"~ paste0("C", C + 5, "H", C*2 - (2*db) +  9, "O8P"),
    class == "PC"  ~ paste0("C", C + 8, "H", C*2 - (2*db) + 16, "NO8P"),
    class == "PE"  ~ paste0("C", C + 5, "H", C*2 - (2*db) + 10, "NO8P"),
    class == "PG"  ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "O10P"),
    class == "PI"  ~ paste0("C", C + 9, "H", C*2 - (2*db) + 15, "O13P"),
    class == "PS"  ~ paste0("C", C + 6, "H", C*2 - (2*db) + 10, "NO10P"),
    
    class == "SM"       ~ paste0("C", C + 5, "H", C*2 - (2*db) + 13,"N2O6P"),
    class == "Cer"      ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO3"),
    class == "Cer;O3"   ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO4"),
    class == "Cer;O4"   ~ paste0("C", C,     "H", C*2 - (2*db) +  1, "NO5"),
    class == "HexCer"   ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO8"),
    class == "HexCer;O3"~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO9"),
    class == "HexCer;O4"~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO10"),
    class == "LactCer"  ~ paste0("C", C +12, "H", C*2 - (2*db) + 21, "NO13"),
    
    class == "MGDG" ~ paste0("C", C +  9, "H", C*2 - (2*db) + 14, "O10"),
    class == "DGDG" ~ paste0("C", C + 15, "H", C*2 - (2*db) + 24, "O15"),
    class == "SQDG" ~ paste0("C", C +  9, "H", C*2 - (2*db) + 14, "O12S"),
    class == "DGTS" ~ paste0("C", C + 10, "H", C*2 - (2*db) + 17, "NO7"),
    
    class == "DGGA" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 12, "O11"),
    
    class == "AI"  ~ paste0("C", 5*M +  9, "H", 10*M - (2*db) + 19, "O11P"),
    class == "AGI" ~ paste0("C", 5*M + 15, "H", 10*M - (2*db) + 29, "O16P"),
    class == "ARC" ~ paste0("C", 5*M +  3, "H", 10*M - (2*db) +  8, "O3"),

    class == "DG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  4, "O5"),
    
    # 3 FA chain -----------------------------------------------------------
    class == "acMGDG"~paste0("C", C +  9, "H", C*2 - (2*db) + 12, "O11"),
    
    class == "TG"   ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O6"),
    class == "TG;O" ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O7"),
    class == "TG;O2" ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O8"),
    class == "TG;O3" ~ paste0("C", C +  3, "H", C*2 - (2*db) +  2, "O9")
  )
}

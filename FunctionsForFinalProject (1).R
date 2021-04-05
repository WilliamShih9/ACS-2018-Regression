

convertList <- function(x, y = c()){
  if ("PUMA" %in% y){
    x$PUMA = convertPUMA(x$PUMA)
    x$PUMA = factor(x$PUMA)
  }
  if ("OCCP" %in% y){
    x$OCCP = convertOCCP(x$OCCP)
    x$OCCP = factor(x$OCCP)
  }
  if ("FOD1P" %in% y){
    x$FOD1P = convertFOD1P(x$FOD1P)
    x$FOD1P = factor(x$FOD1P)
  }
  if ("NAICSP" %in% y){
    x$NAICSP = convertNAICSP(x$NAICSP)
  }
  if ("AGEP" %in% y){
    x$AGEP = convertAGEP(x$AGEP)
  }
  if ("ENG" %in% y){
    x$ENG = convertENG(x$ENG)
  }
  if ("HISP" %in% y){
    x$HISP = convertHISP(x$HISP)
  }
  if ("JWAP" %in% y){
    x$JWAP = convertJWAP(x$JWAP)
  }
  if ("JWMNP" %in% y){
    x$JWMNP = convertJWMNP(x$JWMNP)
  }
  if ("JWTR" %in% y){
    x$JWTR = convertJWTR(x$JWTR)
  }
  if ("MAR" %in% y){
    x$MAR = convertMAR(x$MAR)
  }
  if ("MARHT" %in% y){
    x$MARHT = convertMARHT(x$MARHT)
  }
  if ("QTRBIR" %in% y){
    x$QTRBIR = convertQTRBIR(x$QTRBIR)
  }
  if ("RAC3P" %in% y){
    x$RAC3P = convertRAC3P(x$RAC3P)
  }
  if ("SCHL" %in% y){
    x$SCHL = convertSCHL(x$SCHL)
  }
  if ("SEX" %in% y){
    x$SEX = convertSEX(x$SEX)
  }
  if ("WKHP" %in% y){
    x$WKHP = convertWKHP(x$WKHP)
    x$WKHP = factor(x$WKHP)
  }
  if ("WKW" %in% y){
    x$WKW = convertWKW(x$WKW)
    x$WKW = factor(x$WKW)
  }
  if ("WRK" %in% y){
    x$WRK = convertWRK(x$WRK)
  }
  return(x)
}


regression <- function(x, TheFormula, w = "PWGTP"){
  form = formula(TheFormula)
  gc()
  S = lm(form, weights = get(w), data = x)
  gc()
  AIC.model = AIC(S)
  BIC.model = BIC(S)
  Rsquared = summary(S)[10]
  degrees = summary(S)[[11]][2]
  L = as.data.frame(cbind(AIC.model, BIC.model, Rsquared, degrees))
  colnames(L) = c("AIC","BIC","AdjustedRsquared","df")
  return(L)
}

GetSpeedLM <- function(x, TheFormula, w= "PWGTP"){
    n = length(x[,1])
    one = x[1:400000,]
    two = x[400001:800000,]
    three = x[800001:n,]
    rm(x)
    gc()
    library(speedglm)
    form = formula(TheFormula)
    S = eval(bquote(speedlm(.(form), weights = one[[w]], data = one)))
    gc()
    S = update(S, sparse = FALSE, weights = two[[w]], data = two)
    gc()
    S = update(S, sparse = FALSE, weights = three[[w]], data = three)
    return(S)
}

regressionBig <- function(x, TheFormula, w = "PWGTP"){
    n = length(x[,1])
    one = x[1:400000,]
    two = x[400001:800000,]
    three = x[800001:n,]

    gc()
    library(biglm)
    form = formula(TheFormula)
    S = biglm(form, weights=~get(w), data = one)
    gc()
    S = update(S, two)
    gc()
    S = update(S, three)
    gc()

    AIC.model = n*log(deviance(S)/n) + 2 * length(coef(S))
    BIC.model = n*log(deviance(S)/n) + log(n) * length(coef(S))
    degrees = length(S[[6]])
    Rsquared = summary(S)[[4]]
    ARsquared = 1 - (1 - Rsquared) * ((n-1)/(n - degrees))
    L  = as.data.frame(cbind(AIC.model, BIC.model, ARsquared, degrees))
    colnames(L) = c("AIC","BIC","AdjustedRsquared","df")
    return(L)
}

regressionSpeed <- function(x, TheFormula, w = "PWGTP"){
    n = length(x[,1])
    one = x[1:400000,]
    two = x[400001:800000,]
    three = x[800001:n,]
    gc()
    library(speedglm)
    form = formula(TheFormula)
    S = eval(bquote(speedlm(.(form), weights = one[[w]], data = one)))
    gc()
    S = update(S, sparse = FALSE, weights = two[[w]], data = two)
    gc()
    S = update(S, sparse = FALSE, weights = three[[w]], data = three)
    AIC.model = AIC(S)
    BIC.model = BIC(S)
    degrees = summary(S)[[13]][2]
    Rsquared = summary(S)[[11]]
    L  = as.data.frame(cbind(AIC.model, BIC.model, Rsquared, degrees))
    colnames(L) = c("AIC","BIC","AdjustedRsquared","df")
    return(L)
}


regressionGroup <- function(x, Formula, groups = c(20), variable = "OCCP", type = "PERNP", w = "PWGTP"){
    TheFormula = paste0(Formula,"+","1")
    Frame = regressionBig(x, TheFormula, w)
    x[paste0(variable,0)] = get(paste0("convert",variable))(x[[variable]])
    TheFormula = paste0(Formula,"+",variable,"0")
    Frame[2,] = regressionBig(x, TheFormula, w)
    Number = c(0, length(unique(x[[paste0(variable,0)]]))-1)
    name = c("Null","By Major Group")
    for (i in groups){
        gc()
        name = c(name, paste0(i," equal groups"))
        x[paste0(variable,i+1)] = 
            get(paste0("convert",variable,"quantiles"))(x, type, w, seq(0, 1, length = i+1))
        TheFormula = paste0(Formula,"+",variable,i+1)
        Frame[nrow(Frame)+1,] = regressionBig(x, TheFormula, w)
        Number = c(Number, length(unique(x[[paste0(variable,i+1)]]))-1)
    }
    Frame = cbind(Frame, Number)
    
    gc()
    Frame[,1] = formatC(as.numeric(Frame[,1]), format = "f", digits = 0)
    Frame[,2] = formatC(as.numeric(Frame[,2]), format = "f", digits = 0)
    Frame[,3] = formatC(as.numeric(Frame[,3]), format = "f", digits = 5)
    Frame =
        mutate(Frame,
            AIC = cell_spec(AIC, "latex", bold = (as.numeric(AIC) == min(as.numeric(Frame$AIC)))),
            BIC = cell_spec(BIC, "latex", bold = (as.numeric(BIC) == min(as.numeric(Frame$BIC)))),
            AdjustedRsquared = cell_spec(AdjustedRsquared, "latex",
            bold = (as.numeric(AdjustedRsquared) == max(as.numeric(Frame$AdjustedRsquared))))
        ) 
    Formula = gsub("~"," tilde ",Formula)
    rownames(Frame) = name
    colnames(Frame) = c("AIC", "BIC", "Adjusted R-squared","df","Additional df")
    capt = paste0("AIC, BIC, and Adjusted R-squared of Linear Regression Model for 
    ",Formula,"+",variable," with Different Methods of Splitting Up ",variable)

    PrintTable = kable(Frame, 
             caption = capt, 
             booktabs = TRUE, escape = FALSE, linesep = "") %>% 
             kable_styling(latex_options = c("HOLD_position"), full_width = FALSE) 
    return(PrintTable)
}

regressionList <- function(x, Formulas, title, w = "PWGTP"){
    Frame = regressionBig(x, Formulas[1], w)
    if (length(Formulas) > 1){
        for (i in c(2:length(Formulas))){
            gc()
            Frame[nrow(Frame)+1,] = regressionBig(x, Formulas[i], w)
        }
    }

    Frame[,1] = formatC(as.numeric(Frame[,1]), format = "f", digits = 0)
    Frame[,2] = formatC(as.numeric(Frame[,2]), format = "f", digits = 0)
    Frame[,3] = formatC(as.numeric(Frame[,3]), format = "f", digits = 5)
    Frame =
        mutate(Frame,
            AIC = cell_spec(AIC, "latex", bold = (as.numeric(AIC) == min(as.numeric(Frame$AIC)))),
            BIC = cell_spec(BIC, "latex", bold = (as.numeric(BIC) == min(as.numeric(Frame$BIC)))),
            AdjustedRsquared = cell_spec(AdjustedRsquared, "latex",
            bold = (as.numeric(AdjustedRsquared) == max(as.numeric(Frame$AdjustedRsquared))))
        )    
    colnames(Frame) = c("AIC", "BIC", "Adjusted R-squared","df")  
    rownames(Frame) = sapply(Formulas, function(x) paste0(substr(x, nchar(title)+1, nchar(x))))
    rownames(Frame) = sapply(rownames(Frame), function(x) ifelse(x == "", "Null", x))    
    title = gsub("~"," tilde ",title)
    capt = paste0("AIC, BIC, and Adjusted R-squared of Regression Models in addition to 
    ", title)
    PrintTable = kable(Frame, 
             caption = capt, 
             booktabs = TRUE, escape = FALSE, linesep = "") %>% 
             kable_styling(latex_options = c("HOLD_position"), full_width = FALSE) 
    return(PrintTable)
}

GetFormulaVariations <- function(x, y){
    List = c(1:8)
    List[1] = x
    List[2] = paste0(x,"+",y,"F")
    List[3] = paste0(x,"+","factor(",y,")")
    List[4] = paste0(x,"+",y)
    List[5] = paste0(x,"+",y,"+I(",y,"**2)")
    List[6] = paste0(x,"+",y,"+I(",y,"**2)","+I(",y,"**3)")
    List[7] = paste0(x,"+","sqrt(",y,")+",y)
    List[8] = paste0(x,"+","sqrt(",y,")+",y,"+I(",y,"**2)")
    return(List)
}
    

GetVariables <- function(x = c()){
  Vars = c("PERNP","NAICSP","OCCP","JWAP","WKW","WRK","MAR","MARHT",
           "ENG","JWTR","QTRBIR","RAC3P","SEX","WKHP","SCHL",
           "HISP","JWMNP","AGEP","PUMA","FOD1P","PWGTP")
  Labels = c("Earnings","Industry","Occupation","Time Arrival to Work","Weeks Worked Past Year",
             "Worked Last Week","Marital Status","Times Married","English Fluency",
             "Mode of Transport to Work","Quarter of Birth","Race","Sex",
             "Usual Hours Worked per Week","Educational Attainment",
             "Hispanic origin","Travel Time to Work","Age","Area",
             "Field of Degree","Final Weights")
  A = data.frame(cbind(Vars, Labels))
  if (length(x) > 0){
    A = filter(A, Vars %in% x)
  }
  colnames(A) = c("Variables","Definitions")
  return(A)
}

GetWeights <- function(x = c()){
  Vars = paste0("PWGTP",1:80)
  return(Vars)
}

convertNAICSP <- function(x){
  x[x == ""] = "0"
  S = unique(x)
  None = S[startsWith(S, "0") | startsWith(S,"999920")]
  AGR = S[startsWith(S, "1")]
  EXT = S[startsWith(S, "21")]
  UTL = S[startsWith(S, "22")]
  CON = S[startsWith(S, "23")]
  MFG = S[startsWith(S, "3")]
  WHL = S[startsWith(S, "42")]
  RET = S[startsWith(S, "44") | startsWith(S,"45") | startsWith(S,"4MS")]
  TRN = S[startsWith(S, "48") | startsWith(S,"49")]
  INF = S[startsWith(S, "51")]
  FIN = S[startsWith(S, "52") | startsWith(S,"53")]
  PRF = S[startsWith(S, "54") | startsWith(S, "55") | startsWith(S, "56")]
  EDU = S[startsWith(S, "61")]
  MED = S[startsWith(S, "621") | startsWith(S, "622") | startsWith(S, "623")]
  SCA = S[startsWith(S, "624")]
  ENT = S[startsWith(S, "7")]
  SRV = S[startsWith(S, "8")]
  ADM = S[startsWith(S, "921") | startsWith(S, "923") | startsWith(S, "92M") | 
            startsWith(S, "928P")]
  MIL = S[startsWith(S, "9281")]
  x[x %in% None] = "None"
  x[x %in% AGR] = "AGR"
  x[x %in% EXT] = "EXT"
  x[x %in% UTL] = "UTL"
  x[x %in% CON] = "CON"
  x[x %in% MFG] = "MFG"
  x[x %in% WHL] = "WHL"
  x[x %in% RET] = "RET"
  x[x %in% TRN] = "TRN"
  x[x %in% INF] = "INF"
  x[x %in% FIN] = "FIN"
  x[x %in% PRF] = "PRF"
  x[x %in% EDU] = "EDU"
  x[x %in% MED] = "MED"
  x[x %in% SCA] = "SCA"
  x[x %in% ENT] = "ENT"
  x[x %in% SRV] = "SRV"
  x[x %in% ADM] = "ADM"
  x[x %in% MIL] = "MIL"
  gc()
  x = as.factor(x)
  return(x)
}

PrintNAICSPnames <- function(x){
    industries = c("None","Agriculture, Forestry, Fishing and Hunting",
    "Mining, Quarrying, and Oil and Gas Extraction","Utilities","Construction","Manufacturing",
    "Wholesale Trade","Retail Trade","Transportation and Warehousing",
    "Information","Finance, Insurance, Real Estate, Rental, Leasing",
    "Professional, Scientific, and Technical Services+Management of Companies",
    "Educational Services","Health Care","Social Assistance",
    "Arts, Entertainment, Recreation, Accommodation, Food Services",
    "Other Services","Public Administration","Military")
    name = c("None","AGR","EXT","UTL","CON","MFG","WHL","RET","TRN","INF","FIN","PRF","EDU","MED","SCA","ENT","SRV",
    "ADM","MIL")
    Frame = as.data.frame(cbind(name, industries))
    colnames(Frame) = c("Code","2017 NAICS Sector")
    PrintTable = kable(Frame, 
             caption = "Converting 2018 ACS Industry Code to 2017 NAICS Sectors", 
             booktabs = TRUE, linesep = "") %>% 
             column_spec(1, width = "5em") %>%
             kable_styling(latex_options = c("HOLD_position"), full_width = FALSE)
    return(PrintTable)
}

PrintOCCPnames <- function(x){   
    occupations = c("None","Management","Business","Financial Operations","Computer and Mathematical",
    "Architecture and Engineering","Life, Physical, and Social Science","Community and Social Service",
    "Legal","Education, Training, and Library","Arts, Design, Entertainment, Sports,and Media",
    "Healthcare Practitioners and Technical","Healthcare Support","Protective Service",
    "Food Preparation and Serving Related","Building and Grounds Cleaning and Maintenance",
    "Personal Care and Service","Sales and Related","Office and Administrative Support",
    "Farming, Fishing,and Forestry","Construction and Extraction","Installation, Maintenance, and Repair",
    "Production","Transportation and Material Moving","Military")
    name = c("None","MGR","BUS","FIN","CMM","ENG","SCI","CMS","LGL","EDU","ENT","MED","HLS",
           "PRT","EAT","CLN","PRS","SAL","OFF","FFF","CON","RPR","PRD","TRN","MIL")
    Frame = as.data.frame(cbind(name, occupations))
    colnames(Frame) = c("Code","Major Occupational Group")
    PrintTable = kable(Frame, 
             caption = "Converting 2018 ACS Occupation Code to SOC Major Occupational Groups", 
             booktabs = TRUE, linesep = "") %>% 
             column_spec(1, width = "5em") %>%
             kable_styling(latex_options = c("HOLD_position"), full_width = FALSE)
    return(PrintTable)
}

convertOCCP <- function(x){
  x[x %in% 9920] = 0
  occupations = c(-1, 1, 440,750,960,1240,1560,1980,2060,2180,2555,2920,3550,3655,3960,4160,
                  4255,4655,4965,5940,6130,6950,7640,8990,9760,9830)
  name = c("None","MGR","BUS","FIN","CMM","ENG","SCI","CMS","LGL","EDU","ENT","MED","HLS",
           "PRT","EAT","CLN","PRS","SAL","OFF","FFF","CON","RPR","PRD","TRN","MIL")
  x = cut(x, breaks = occupations)
  levels(x) = name
  return(x)
}


convertNAICSPquantiles <- function(x, type = "PERNP", w = "PWGTP", quant = seq(0, 1, 0.05)){
  x = select(x, w, "NAICSP", type)
  y = group_split(x, NAICSP)
  k = group_keys(x, NAICSP)
  m = sapply(y, weightedmean, type, w)
  p = sapply(y, weight, w)
  combined = cbind(k, m, p)
  combined = arrange(combined, m)
  combined[,3] = cumsum(combined[,3])/sum(combined[,3])
  combined$p = cut(combined$p, breaks = quant)
  Split = group_split(combined, p)
  Groups = sapply(Split, function(x) x[,1]) 
  l = as.character(unique(combined$p))
  l = strsplit(substring(l, first = 2, last = nchar(l)-1), ",")
  l = sapply(l, function(x) paste0(as.numeric(x[1])*100,"%-",as.numeric(x[2])*100, "%"))
  x = x$NAICSP 
  for (i in c(1:length(Groups))){
    x[x %in% Groups[[i]]] = l[i]
  }
  return(as.factor(x))
}

convertOCCPquantiles <- function(x, type = "PERNP", w = "PWGTP", quant = seq(0, 1, 0.05)){
  x = select(x, w, "OCCP", type)
  y = group_split(x, OCCP)
  k = group_keys(x, OCCP)
  m = sapply(y, weightedmean, type, w)
  p = sapply(y, weight, w)
  combined = cbind(k, m, p)
  combined = arrange(combined, m)
  combined[,3] = cumsum(combined[,3])/sum(combined[,3])
  combined$p = cut(combined$p, breaks = quant)
  Split = group_split(combined, p)
  Groups = sapply(Split, function(x) x[,1]) 
  l = as.character(unique(combined$p))
  l = strsplit(substring(l, first = 2, last = nchar(l)-1), ",")
  l = sapply(l, function(x) paste0(as.numeric(x[1])*100,"%-",as.numeric(x[2])*100, "%"))
  x = x$OCCP 
  for (i in c(1:length(Groups))){
    x[x %in% Groups[[i]]] = l[i]
  }
  return(as.factor(x))
}


convertJWAP <- function(x){
  times = c(-1,0,45,57,63,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,
            120,123,129,141,153,165,213, 285)
  name = c("Home","12AM-4AM","4AM-5AM","5AM-5:30AM","5:30AM-6AM","6:00AM-6:15AM","6:15AM-6:30AM","6:30AM-6:45AM",
           "6:45AM-7AM","7AM-7:15AM","7:15AM-7:30AM","7:30AM-7:45AM","7:45AM-8AM",
           "8AM-8:15AM","8:15AM-8:30AM","8:30AM-8:45AM","8:45AM-9AM","9AM-9:15AM",
           "9:15AM-9:30AM","9:30AM-9:45AM","9:45AM-10AM","10AM-10:15AM","10:15AM-10:30AM","10:30AM-11AM",
           "11AM-12PM","12PM-1PM","1PM-2PM","2PM-6PM","6PM-11:59PM")
  x = cut(x, breaks = times)
  levels(x) = name
  return(x)
}

convertWKW <- function(x){
    weeks = c(-1,0,1,2,3,4,5,6)
    name = c("0wks","50-52wks","48-49wks","40-47wks","27-39wks","14-26wks","01-13wks")
    x = cut(x, breaks = weeks)
    levels(x) = name
    return(x)
}

convertWRK <- function(x){
    work = c(-1,0,1,2)
    name = c("Not reported","Worked last week","Did not work last week")
    x = cut(x, breaks = work)
    levels(x) = name
    return(x)
}


convertMAR <- function(x){
    marriagestatus = c(0,1,2,3,4,5)
    name = c("Married","Widowed","Divorced","Separated","Never Married")
    x = cut(x, breaks = marriagestatus)
    levels(x) = name
    return(x)
}

convertMARHT <- function(x){
    marriagetimes = c(-1,0,1,2,3)
    name = c("Never Married","Married Once","Married Twice","Married 3+ Times")
    x = cut(x, breaks = marriagetimes)
    levels(x) = name
    return(x)
}

convertENG <- function(x){
    englishlit = c(-1,0,1,2,3,4)
    name = c("Only English","Very well","Well","Not well","Not at all")
    x = cut(x, breaks = englishlit)
    levels(x) = name
    return(x)
}

convertJWTR <- function(x){
    commutes = c(-1,0,1,2,3,4,5,6,7,8,9,10,12)
    name = c("Not at work","Car/truck/van","Bus","Streetcar","Subway","Railroad","Ferryboat","Taxi","Motorcycle",
        "Bicycle","Walked","Worked at home/Other")
    x = cut(x, breaks = commutes)
    levels(x) = name
    return(x)
}

convertQTRBIR <- function(x){
    quarter = c(0,1,2,3,4)
    name = c("Jan-Mar","Apr-June","July-Sept","Oct-Dec")
    x = cut(x, breaks = quarter)
    levels(x) = name
    return(x)
}

    
    

convertRAC3P <- function(x){
  races = c(0,1,2,3,4,5,6,14,15,100)
  name = c("White Alone","Black Alone","Native American Alone","Indian alone","Chinese alone",
           "Filipino alone","Other Asian or Pacific Islander","Some Other Race alone","Mixed Race")
  x = cut(x, breaks = races)
  levels(x) = name
  return(x)
}


convertWKHP <- function(x){
  hours = c(-1,0,5,10,15,20,25,30,35,39,40,45,50,55,60,65,70,99)
  name = c("0hrs","01-5hrs","06-10hrs","11-15hrs","16-20hrs","21-25hrs","26-30hrs",
           "31-35hrs","36-39hrs","40hrs","41-45hrs","46-50hrs","51-55hrs","56-60hrs","61-65hrs",
           "66-70hrs",">70hrs")
  x = cut(x, breaks = hours)
  levels(x) = name
  return(x)
}

convertSEX <- function(x){
    sexes = c(0,1,2)
    name = c("Male","Female")
    x = cut(x, breaks = sexes)
    levels(x) = name
    return(x)
}

#Convert educational attainment to factor
convertSCHL <- function(x){
  degrees = c(0,15,17,19,20,21,22,23,24)
  name = c("Less than HS diploma","HS or equivalent","Some college","Associate's","Bachelor's",
        "Master's","Professional","Doctorate")
  x = cut(x, breaks = degrees)
  levels(x) = name
  return(x)
}

#Convert Hispanic origin to factor
convertHISP <- function(x){
  hispanicorigin = c(0,1,2,24)
  name = c("Not Hispanic","Mexican","Hispanic,not Mexican")
  x = cut(x, breaks = hispanicorigin)
  levels(x) = name
  return(x)
}


#Convert Commute Time to factor
convertJWMNP <- function(x){
  commutes = c(-1,1,5,10,15,20,25,30,35,40,45,50,55,60,80,100,999)
  name = c("00-01min","02-5min","06-10min","11-15min","16-20min","21-25min","26-30min","31-35min",
           "36-40min","41-45min","46-50min","51-55min","56-60min","61-80min",
           "81-100min",">100min")
  x = cut(x, breaks = commutes)
  levels(x) = name
  return(x)
}

#Convert Age to factor
convertAGEP <- function(x){
  ages = c(24,25,26,27,28,29,30,31,33,35,37,39,44,54,59,64)
  name = c("25","26","27","28","29","30","31","32-33","34-35","36-37","38-39","40-44","45-54","55-59","60-64")
  x = cut(x, breaks = ages)
  levels(x) = name
  return(x)
}



#Convert PUMA to quantiles of mean income of PUMA
convertPUMAquantiles <- function(x, type = "PERNP", w = "PWGTP", quant = seq(0, 1, 0.05)){
  x = select(x, "PUMA", w, type)
  Temp = x$PUMA
  Split = group_split(x, PUMA)
  Med = cbind(group_keys(x, PUMA), sapply(Split, weightedmean, type))
  names(Med) = c("PUMA",type)
  Quants = quantile(Med[,2], quant)
  Quants = Quants[!duplicated(Quants)]
  Med[,2] = cut(Med[,2], breaks = Quants)
  Split2 = group_split(Med, get(type))
  PUMAsplit = sapply(Split2, function(x) return(x[,1]))
  QName = names(Quants)
  Name = c()
  for (i in c(1:length(Quants)-1)){
    Name[i] = paste0(QName[i],"-",QName[i+1])
  }
  if (length(PUMAsplit) > length(Name)){
    PUMAsplit[[1]] = c(PUMAsplit[[1]], PUMAsplit[[length(PUMAsplit)]])
  }
  for (i in c(1:length(PUMAsplit))){
    Temp[Temp %in% PUMAsplit[[i]]] = Name[i]
  }
  Temp = as.factor(Temp)
  return(Temp)
}


#Convert PUMA to arbitrarily decided regions where large city within each state are grouped
convertPUMA <- function(x){
  #Alabama 1
  
  JeffersonAL = c(101301:101305)
  Mobile = c(102701:102703) 
  x[x %in% JeffersonAL] = "AL,Jefferson"
  x[x %in% Mobile] = "AL,Mobile"
  x[x %in% c(100000:102703)] = "Other"
  #Alaska 2
  x[x %in% c(200101:200400)] = "Other"
  #Arizona 3
  Maricopa = c(400100:400134)
  Pima = c(400201:400209)
  x[x %in% Maricopa] = "AZ,Maricopa"
  x[x %in% Pima] = "AZ,Pima"
  x[x %in% c(400135:400900)] = "Other"
  gc()
  #Arkansas 5
  x[x %in% c(500100:502000)] = "Other"
  #California 6
  Alameda = c(600101:600110)
  ContraCosta = c(601301:601309)
  Fresno = c(601901:601907)
  Kern = c(602901:602905)
  LosAngeles = c(603701:603769)
  OrangeCA = c(605901:605918)
  Riverside = c(606501:606515)
  Sacramento = c(606701:606712)
  SanBernardino = c(607101:607115)
  SanDiego = c(607301:607322)
  SF = c(607501:607507)
  SanMateo = c(608101:608106)
  SantaClara = c(608501:608514)
  Ventura = c(611101:611106)
  x[x %in% Alameda] = "CA,Alameda"
  x[x %in% ContraCosta] = "CA,ContraCosta"
  x[x %in% Fresno] = "CA,Fresno"
  x[x %in% Kern] = "CA,Kern"
  x[x %in% LosAngeles] = "CA,LA"
  x[x %in% OrangeCA] = "CA,Orange"
  x[x %in% Riverside] = "CA,Riverside"
  x[x %in% Sacramento] = "CA,Sacramento"
  x[x %in% SanBernardino] = "CA,SanBernadino"
  x[x %in% SanDiego] = "CA,SD"
  x[x %in% SF] = "CA,SF"
  x[x %in% SanMateo] = "CA,SM"
  x[x %in% SantaClara] = "CA,SC"
  x[x %in% Ventura] = "CA,Ventura"
  x[x %in% c(600300:611300)] = "Other"
  #Colorado 8
  Denver = c(800812:800816)
  ColoradoSprings = c(804101:804106)
  DenverSuburb = c(800801:800811,800817:800824)
  x[x %in% Denver] = "CO,Denver"
  x[x %in% ColoradoSprings] = "CO,ElPaso"
  x[x %in% DenverSuburb] = "CO,DenverS"
  x[x %in% c(800100:804106)] = "Other"
  #Connecticut 9
  NewHaven = c(900901:900906)
  Hartford = c(900300:900306)
  x[x %in% NewHaven] = "CT,New"
  x[x %in% Hartford] = "CT,Hart"
  x[x %in% c(900100:901500)] = "Other"
  #Delaware 10
  x[x %in% c(1000101:1000300)] = "DE,Delaware"
  #District of Columbia 11
  x[x %in% c(1100101:1100105)] = "Urban"
  #Florida 12
  Broward = c(1201101:1201114)
  Duval = c(1203101:1203107)
  Hillsborough = c(1205701:1205708)
  MiamiDade = c(1208601:1206824)
  OrangeFL = c(1209501:1209510)
  PalmBeach = c(1209901:1209911)
  Pinellas = c(1210301:1210308)
  x[x %in% Broward] = "FL,Broward"
  x[x %in% Duval] = "FL,Duval"
  x[x %in% Hillsborough] = "FL,Hillsbourough"
  x[x %in% OrangeFL] = "FL,Orange"
  x[x %in% PalmBeach] = "FL,PalmBeach"
  x[x %in% Pinellas] = "FL,Pinellas"
  x[x %in% c(1200101:1212704)] = "Other"
  #Georgia 13
  CAtlanta = c(1301001:1301008)
  x[x %in% CAtlanta] = "GA,Atlanta"
  x[x %in% c(1300100:1306002)] = "Other"
  #Hawaii 15
  x[x %in% c(1500100:1500308)] = "HI,Hawaii"
  #Idaho 16
  x[x %in% c(1600100:1601300)] = "Other"
  #Illinois 17
  Cook = c(1703401:1703422)
  Chicago = c(1703501:1703532)
  x[x %in% Cook] = "Urban"
  x[x %in% Chicago] = "Urban"
  x[x %in% c(1700104:1703700)] = "Other"
  #Indiana 18
  Marion = c(1802301:1802307)
  x[x %in% Marion] = "IN,Marion"
  x[x %in% c(1800101:1803600)] = "Other"
  #Iowa 19
  x[x %in% c(1900100:1902300)] = "IA,Iowa"
  #Kansas 20
  x[x %in% c(2000100:2001500)] = "KS,Kansas"
  #Kentucky 21
  Jefferson = c(2101701:2101706)
  x[x %in% Jefferson] = "KY,Jefferson"
  x[x %in% c(2100100:2102800)] = "Other"
  #Louisiana 22
  NewOrleans = c(2202300:2202302,2202400:2202402)
  x[x %in% NewOrleans] = "LA,NewOrleans"
  x[x %in% c(2200100:2202500)] = "Other"
  #Maine 23
  x[x %in% c(2300100:2301000)] = "ME,Maine"
  #Maryland 24
  BaltimoreCounty = c(2400501:2400507)
  BaltimoreCity = c(2400801:2400805)
  PrinceGeorge = c(2401101:2401107)
  Montgomery = c(2401001:2401005)
  x[x %in% BaltimoreCounty] = "Baltimore"
  x[x %in% BaltimoreCity] = "BaltimoreCity"
  x[x %in% PrinceGeorge] = "PrinceGeorge"
  x[x %in% Montgomery] = "MD,Montgomery"
  x[x %in% c(2400100:2401600)] = "Other"
  #Massachusetts 25
  Middlesex = c(2500501:2500508)
  Boston = c(2503301:2503306)
  x[x %in% Middlesex] = "Middlesex"
  x[x %in% Boston] = "Boston"
  x[x %in% c(2500100:2504903)] = "Other"
  #Michigan 26
  Oakland = c(2602901:2602908)
  Wayne = c(2603201:2603213)
  x[x %in% Oakland] = "MI,Oakland"
  x[x %in% Wayne] = "MI,Wayne"
  x[x %in% c(2600100:2603300)] = "Other"
  #Minnesota 27
  Hennepin = c(2701401:2701410)
  x[x %in% Hennepin] = "MN,Hennepin"
  x[x %in% c(2700100:2702600)] = "Other"
  #Mississippi 28
  x[x %in% c(2800100:2802100)] = "MS,Missisippi"
  #Missouri 29
  StLouis = c(2901801:2901808)
  x[x %in% StLouis] = "MO,StLouis"
  x[x %in% c(2900100:2902800)] = "Other"
  #Montana 30
  x[x %in% c(3000100,3000200,3000300,3000400,3000500,3000600,3000700)] = "MT,Montana"
  #Nebraska 31
  x[x %in% c(3100100,3100200,3100300,3100400,3100500,3100600,3100701:3100904)] = "NE,Nebraska"
  #Nevada 32
  Clark = c(3200401:3200413)
  Other = c(3200101:3200300)
  x[x %in% Clark] = "NV,Clark"
  x[x %in% Other] = "Other"
  #New Hampshire 33
  x[x %in% c(3300100:3301000)] = "NH,NewHampshire"
  #New Jersey 34
  Bergen = c(3400301:3400308)
  Middlesex = c(3400901:3400907)
  x[x %in% Bergen] = "NJ,Bergen"
  x[x %in% Middlesex] = "NJ,Middlesex"
  x[x %in% c(3400101:3402600)] = "Other"
  #New Mexico 35
  Albuquerque = c(3500801:3500806)
  x[x %in% Albuquerque] = "NM,Albuquerque"
  x[x %in% c(3500100:3501200)] = "Other"
  #New York 36
  Erie = c(3601201:3601207)
  Nassau = c(3603201:3603212)
  Suffolk = c(3603301:3603313)
  Bronx = c(3603701:3603710)
  Manhattan = c(3603801:3603810)
  Brooklyn = c(3604001:3604018)
  Queens = c(3604101:3604114)
  x[x %in% Erie] = "NY,Erie"
  x[x %in% Nassau] = "NassauNY"
  x[x %in% Suffolk] = "SuffolkNY"
  x[x %in% Bronx] = "Bronx"
  x[x %in% Manhattan] = "Manhattan"
  x[x %in% Brooklyn] = "BrooklynNY"
  x[x %in% Queens] = "QueensNY"
  x[x %in% c(3600100:3604114)] = "Other"
  #North Carolina 37
  Wake = c(3701201:3701208)
  Mecklenburg = c(3703101:3703108)
  x[x %in% Wake] = "NC,Wake"
  x[x %in% Mecklenburg] = "NC,Mecklenburg"
  x[x %in% c(3700100:3705400)] = "Other"
  #North Dakota 38
  gc()
  x[x %in% c(3800100,3800200,3800300,3800400,3800500)] = "Other"
  #Ohio 39
  Cuyahoga = c(3900901:3900910)
  Columbus = c(3904101:3904111)
  Hamilton = c(3905501:3905507)
  x[x %in% Cuyahoga] = "OH,Cuya"
  x[x %in% Columbus] = "OH,Colum"
  x[x %in% Hamilton] = "OH,Hamilton"
  x[x %in% c(3900100:3905700)] = "Other"
  #Oklahoma 40
  Oklahoma = c(4001001:4001006)
  x[x %in% Oklahoma] = "OK,OklahomaCity"
  x[x %in% c(4000100:4001601)] = "Other"
  #Oregon 41
  PortlandCity = c(4101301:4101316)
  x[x %in% PortlandCity] = "Urban"
  x[x %in% c(4100100:4101324)] = "Other"
  #Pennsylvania 42
  Allegheny = c(4201801:4201807)
  Montgomery = c(4203101:4203106)
  Philadelphia = c(4203201:4203211)
  x[x %in% Allegheny] = "PA,Allegheny"
  x[x %in% Montgomery] = "PA,Montgomery"
  x[x %in% Philadelphia] = "PA,Philly"
  x[x %in% c(4200101:4204002)] = "Other"
  #Puerto Rico
  x[x %in% c(7200101:7201102)] = "Puerto Rico"
  #Rhode Island 44
  x[x %in% c(4400101:4400400)] = "RI,Rhode Island"
  #South Carolina 45
  x[x %in% c(4500101:4501600)] = "SC,South Carolina"
  #South Dakota 46
  x[x %in% c(4600100,4600200,4600300,4600400,4600500,4600600)] = "Other"
  #Tennessee 47
  Shelby = c(4703201:4703208)
  x[x %in% Shelby] = "TN,Shelby"
  x[x %in% c(4700100:4703208)] = "Other"
  #Texas 48
  Dallas = c(4802301:4802322)
  Tarrant = c(4802501:4802516)
  Houston = c(4804601:4804638)
  Austin = c(4805301:4805309)
  SanAntonio = c(4805901:4805916)
  x[x %in% Dallas] = "TX,Dalas"
  x[x %in% Tarrant] = "TX,Tarrant"
  x[x %in% Houston] = "TX,HOuston"
  x[x %in% Austin] = "TX,AUstin"
  x[x %in% SanAntonio] = "TX,SA"
  x[x %in% c(4800100:4806900)] = "Other"
  #Utah 49
  SaltLake = c(4935001:4935009)
  x[x %in% SaltLake] = "UT,SaltLake"
  x[x %in% c(4903001:4957002)] = "Other"
  #Vermont 50
  x[x %in% c(5000100,5000200,5000300,5000400)] = "Other"
  #Virginia 51
  Fairfax = c(5159301:5159309)
  x[x %in% Fairfax] = "VA,Fairfax"
  x[x %in% c(5101301:5159309)] = "Other"
  #Washington 52
  Pierce = c(5311501:5311507)
  King = c(5311601:5311616)
  x[x %in% Pierce] = "WA,Pierce"
  x[x %in% King] = "WA,King"
  x[x %in% c(5310100:5311900)] = "Other"
  #West Virginia 54
  x[x %in% c(5400100:5401300)] = "WV,West Virginia"
  #Wisconsin 55
  Milwaukee = c(5540101:5541005)
  x[x %in% Milwaukee] = "WI,Milwaukee"
  x[x %in% c(5500100:5570301)] = "Other"
  #Wyoming 56
  x[x %in% c(5600100:5600500)] = "Other"
  x = as.factor(x)
  return(x)
}

#Convert field of degree to factors
convertFOD1P <- function(x){
  x[x %in% 0] = "NoDegree"
  x[x %in% c(1100:1199)] = "Agriculture"
  x[x %in% c(1301:1303)] = "Environment"
  x[x %in% c(2100:2107)] = "Computer"
  x[x %in% c(2300:2399)] = "Education"
  x[x %in% c(1401,2400:2599)] = "Engineering"
  x[x %in% c(2600:2603,3301:3302)] = "Language"
  x[x %in% c(3600:3699)] = "Biology"
  x[x %in% c(3700:3702)] = "Math"
  x[x %in% c(4000:4007)] = "Multiple"
  x[x %in% c(5000:5098)] = "Physical Sciences"
  x[x %in% c(5200:5299)] = "Psychology"
  x[x %in% c(5500:5599)] = "Social Science"
  x[x %in% c(6000:6099)] = "Art"
  x[x %in% c(6100:6199)] = "Health"
  x[x %in% c(6200:6299)] = "Business"
  x[x %in% c(1:9999)] = "Other"
  x = as.factor(x)
  return(x)
}


convertFOD1Pquantiles <- function(x, type = "PERNP", w = "PWGTP", quant = seq(0, 1, 0.05)){
  x = select(x, w, "FOD1P", type)
  y = group_split(x, FOD1P)
  k = group_keys(x, FOD1P)
  m = sapply(y, weightedmean, type, w)
  p = sapply(y, weight, w)
  combined = cbind(k, m, p)
  combined = arrange(combined, m)
  combined[,3] = cumsum(combined[,3])/sum(combined[,3])
  combined$p = cut(combined$p, breaks = quant)
  Split = group_split(combined, p)
  Groups = sapply(Split, function(x) x[,1]) 
  l = as.character(unique(combined$p))
  l = strsplit(substring(l, first = 2, last = nchar(l)-1), ",")
  l = sapply(l, function(x) paste0(as.numeric(x[1])*100,"%-",as.numeric(x[2])*100, "%"))
  x = x$FOD1P 
  for (i in c(1:length(Groups))){
    x[x %in% Groups[[i]]] = l[i]
  }
  x = as.factor(x)
  return(x)
}


#Calculates the estimates and the standard errors of the regression coefficients
SEregression <- function(x, TheFormula){
  form = formula(TheFormula)
  replicates = paste0("PWGTP",1:80)
  regression = c()
  fin = summary(lm(form, weights = PWGTP, data = x))[5:11]
  final = fin[[1]][,1]
  l = vector(mode = "list", length = 81)
  l[[1]] = fin
  gc()
  df = summary(lm(form, weights = PWGTP, data = x))[[8]][2]
  for (i in c(1:80)){
    weigh = x[replicates[i]]
    weigh[weigh < 0] = 0
    intermed = summary(lm(form, weights = weigh[[1]], data = x))[5:11]
    gc()
    l[[i+1]] = intermed
    regression = rbind(regression, intermed[[1]][,1])
    gc()
  }
  s = 0
  for (i in c(1:80)){
    s = s + (regression[i,] - final)^2
  }
  s = (s * 4/80)^0.5
  s = cbind(final, s, final/s, 2*pt(final/s, df, lower = FALSE))
  colnames(s) = c("Estimate","Standard Error","t value","Pr(>|t|)")
  return(list(s,l))
}

#These means are gotten from ACS 2007, when the last time that detailed hours worked was available
#The weighted mean of each week range (from ACS 2018) was taken and used 
#Gets wage per hour from weeks worked, hours worked per week, and total earnings
WagePerHour <- function(x){
  Weeks = x$WKW
  Weeks[Weeks == 1] = 51.85 #50-52
  Weeks[Weeks == 2] = 48.19 #48-49
  Weeks[Weeks == 3] = 42.38 #40-47
  Weeks[Weeks == 4] = 33.08 #27-39
  Weeks[Weeks == 5] = 21.32 #14-26
  Weeks[Weeks == 6] = 7.49 #1-13
  Hours = Weeks*x$WKHP
  Wage = x$PERNP/Hours
  return(Wage)
}


#Calculates the standard errors of the total population of the sample
SEpopulation <- function(x){
  replicates = paste0("PWGTP",1:80)
  final = sum(x[["PWGTP"]])
  eighty = c()
  for(i in c(1:80)){
    eighty = c(eighty, sum(x[replicates[i]]))
  }
  s = 0
  for(i in c(1:80)){
    s = s + (eighty[i] - final)^2
  }
  s = (s * 4/80)^0.5
  return(s)
}

#Calculates the standard errors of the mean of one of the variables in the sample
SEmean <- function(x, variable){
  factor = variable
  replicates = paste0("PWGTP",1:80)
  final = weightedmean(x, variable)
  eighty = c()
  for (i in c(1:80)){
    eighty = c(eighty, weightedmean(x, variable, replicates[i]))
  }
  s = 0
  for(i in c(1:80)){
    s = s + (eighty[i] - final)^2
  }
  s = (s * 4/80)^0.5
  return(s)
}

#The population estimate of this sample

weight <- function(x, w = "PWGTP"){
  return(sum(x[w]))
}

#The mean of one of the variables in this sample
weightedmean <- function(x, variable, weights = "PWGTP"){
  one = x[weights]
  two = x[variable]
  average = sum(one*two)/sum(one)
  return(average)
}

#The weighted median of the variable not assuming the variable is sorted low to high
weightedmedian <- function(x, variable, weights = "PWGTP"){
  one = x[weights]
  quantile = sum(one)/2
  two = x[variable]
  combined = as.data.frame(cbind(one, two))
  combined = arrange(combined, get(variable))
  fortypercent = 0.4*length(one[[1]])
  comp = sum(combined[1:fortypercent,1])
  comptemp = comp
  for (i in c(fortypercent:length(one[[1]]))){
    comptemp = comp
    comp = comp + combined[i,1]
    if (comp > quantile){
      break
    }
  } 
  offsetpositive = quantile - comptemp
  offsetnegative = comp - quantile
  offsetpositive = offsetpositive/(comp-comptemp)
  offsetnegative = offsetnegative/(comp-comptemp)
  result = combined[i,2]*offsetpositive + combined[i-1,2]*offsetnegative
  return(result)
}

#Calculates the quantiles of the variable of this sample x. No assumption that the variable is sorted low to high
weightedquantiles <- function(x, variable, quantiles = seq(0.1,0.9,0.1), weights = "PWGTP"){
  one = x[weights]
  quantile = quantiles*sum(one)
  two = x[variable]
  combined = as.data.frame(cbind(one, two))
  combined = arrange(combined, get(variable))
  comp = 0
  count = 1
  ret = c()
  for (i in c(1:length(one[[1]]))){
    comptemp = comp
    comp = comp + combined[i,1]
    if (comp > quantile[count]){ 
      offsetpositive = quantiles[count] - comptemp
      offsetnegative = comp - quantiles[count]
      offsetpositive = offsetpositive/(comp-comptemp)
      offsetnegative = offsetnegative/(comp-comptemp)
      result = combined[i,2]*offsetpositive + combined[i-1,2]*offsetnegative
      ret = c(ret, result)
      count = count + 1
    }
    if (count > length(quantiles)){
      break
    }
  } 
  return(ret)
}


GetSummaryStatistics <- function(x, variable, quantiles = FALSE){
    if (quantiles){
        x[variable] = get(paste0("convert",variable,"quantiles"))(x)
    }
    else{
        x = convertList(x, variable)
    }
    Split = group_split(x, get(variable))
    
    Salary = sapply(Split, weightedmean, "PERNP")
    SDSalary = sapply(Split, weightedmedian, "PERNP")
    Hour= sapply(Split, weightedmean, "PERNPHR")
    SDHour = sapply(Split, weightedmedian, "PERNPHR")
    Weight = sapply(Split, weight)
    
    Tab = cbind(group_keys(x, get(variable)), Salary, SDSalary, Hour, SDHour, Weight)
    Tab = arrange(Tab, desc(Salary))
    Tab[,c(2:5)] = exp(Tab[,c(2:5)])
    
    Tab[,6] = prettyNum(Tab[,6], big.mark = ",", digits = 0, scientific = FALSE)
    Tab[,5] = formatC(Tab[,5], big.mark = ",", format = "f", digits = 2)
    Tab[,4] = formatC(Tab[,4], big.mark = ",", format = "f", digits = 2)
    Tab[,3] = prettyNum(Tab[,3], big.mark = ",", digits = 0, scientific = FALSE)
    Tab[,2] = prettyNum(Tab[,2], big.mark = ",", digits = 0, scientific = FALSE)
    
    Name = as.character(GetVariables(variable)[,2])
    
    colnames(Tab) = c(Name,"Log-Average","Median","Log-Average","Median","# of Workers")
    
    capt = paste("Geometric Mean and Median Earnings per Year/Hour split by", Name)   
    
    PrintTable = kable(Tab, "latex",
             caption = capt, 
             booktabs = TRUE, linesep = "") %>% 
             kable_styling(latex_options = c("HOLD_position"), full_width = FALSE) %>%
             column_spec(1, bold = T) %>%
             column_spec(2, border_left = T) %>%
             column_spec(4, border_left = T) %>%
             column_spec(6, border_left = T) %>%
             add_header_above(c(" " = 1, "Earnings Per Year" = 2, "Earnings Per Hour" = 2, " " = 1))
    return(PrintTable)
}  

GetCoefficientPercent <- function(Model, Term, Term2){
    Model$coefficients = Model$coefficients[,1]
    N = Model$xlevels[[Term2]][1]
    One = data.frame(0.00)
    A = as.data.frame(Model$coefficients[grepl(Term, names(Model$coefficients))])
    rownames(One) = N
    colnames(One) = colnames(A)
    A = rbind(One, A) 
    A = tibble::rownames_to_column(A)
    B = strsplit(A[,1], "-", fixed = TRUE)
    B = data.frame(matrix(unlist(B), nrow = length(B), byrow = TRUE))
    A = cbind(A,B)
    A[,4] = as.character(A[,4])
    A[,4] = as.numeric(gsub("%","",A[,4]))
    A = arrange(A, X2)
    B = cbind(A[,4], A[,2])
    colnames(B) = c("Percent",Term)
    B = data.frame(B)
    return(B)
}


ConfidenceInterval <- function(Model, Term, Term2){
    colnames(Model$coefficients) = c("One","Two","Three","Four")
    N = Model$xlevels[[Term2]][1]
    One = data.frame(cbind(0.00,0,Model$coefficients[1,3],
            Model$coefficients[1,4]))
    A = Model$coefficients[grepl(Term, rownames(Model$coefficients)),] 
    rownames(One) = N
    colnames(One) = colnames(A)
    A = rbind(One, A)
    Keep = rownames(A)
    Upper = A[,1] + 2.576*A[,2]
    Lower = A[,1] - 2.576*A[,2]
    A = as.data.frame(cbind(A[,1], Lower, Upper))
    rownames(A) = gsub(Term, "", Keep)
    A = tibble::rownames_to_column(A)
    colnames(A) = c(Term,"Estimate","99% Lower","99% Upper")
    return(A)
}

PrintConfidenceIntervalTable <- function(H, Y, Term, Term2){
    x = cbind(ConfidenceInterval(H, Term, Term2), ConfidenceInterval(Y, Term, Term2))
    library(scales)
    Hour = x[,3:4]
    Year = x[,7:8]
    Name = x[,1]
    Hour = exp(Hour)
    Year = exp(Year)
    Hour[,1] = percent(Hour[,1])
    Hour[,2] = percent(Hour[,2])
    Year[,1] = percent(Year[,1])
    Year[,2] = percent(Year[,2])
    HourText = paste0("(",Hour[,1],", ",Hour[,2],")")
    YearText = paste0("(",Year[,1],", ",Year[,2],")")
    
    Tab = cbind(HourText, YearText)
    Tab[1,] = c(" Baseline ", " Baseline ")
    capt = paste0("99% Confidence Intervals for ",colnames(x)[1], " relative to ", Name[1])
    colnames(Tab) = c("Earnings per Hour","Earnings per Year")
    rownames(Tab) = Name
    myHeader = c(capt = 3)
    names(myHeader) = c(capt)
    A = kable(Tab, "latex", booktabs = TRUE, linesep = "") %>% 
        kable_styling(Tab, latex_options = c("HOLD_position")) %>%
        column_spec(1, bold = TRUE) %>%
        column_spec(2:3, border_left = TRUE) %>%
        add_header_above(header = myHeader, bold = TRUE)
    return(A)
}

PrintSCHLTable <- function(Model){
    A = Model$coefficients[grepl("SCHL", rownames(Model$coefficients)),]
    Thirty = A[1:7,1] + log2(30)*A[8:14,1]
    Forty= A[1:7,1] + log2(40)*A[8:14,1]
    Fifty = A[1:7,1] + log2(50)*A[8:14,1]
    Sixty = A[1:7,1] + log2(60)*A[8:14,1]
    SEthirty = sqrt(A[1:7,2]^2 + log2(30)*A[8:14,2]^2)
    SEforty = sqrt(A[1:7,2]^2 + log2(40)*A[8:14,2]^2)
    SEfifty = sqrt(A[1:7,2]^2 + log2(50)*A[8:14,2]^2)
    SEsixty = sqrt(A[1:7,2]^2 + log2(0)*A[8:14,2]^2)
    Combined = cbind(Thirty,Forty,Fifty,Sixty)
    CombinedSE = cbind(SEthirty,SEforty,SEfifty,SEsixty)
    Lower = exp(Combined - 2.576*CombinedSE)
    Upper = exp(Combined + 2.576*CombinedSE)
    
    Lower = apply(Lower, 2, percent)
    Upper = apply(Upper, 2, percent)
    Text = matrix(paste0("(", Lower, ", ", Upper, ")"), nrow = 7, ncol = 4, byrow = F)
    Text = rbind("Baseline",Text)
    rownames(Text) = c("Less than HS","HS or equivalent","Some college","Associate's",
    "Bachelor's","Master's","Professional","Doctorate")
    colnames(Text) = c("Age 30","Age 40","Age 50","Age 60")
    capt = "99% Confidence Intervals for Educational Attainment relative to less than HS for Hour model"
    myHeader = c(capt = 5)
    names(myHeader) = c(capt)
    A = knitr::kable(Text, "latex", booktabs = TRUE, linesep = "") %>%
        kable_styling(latex_options = c("HOLD_position")) %>%
        column_spec(1, bold = TRUE) %>%
        column_spec(2:5, border_left = TRUE) %>%
        add_header_above(header = myHeader, bold = TRUE)
    return(A)
}

PrintMARTable <- function(Model){
    A = Model$coefficients[grepl("SEX|MAR", rownames(Model$coefficients)),]
    B = Model$cov[grepl("MAR",rownames(Model$coefficients)),
                grepl("SEXFemale",rownames(Model$coefficients))][1:4]
    MaleAGEP = log2(25)*A[6,1]
    FemaleAGEP = log2(25)*A[7,1]
    Female = FemaleAGEP + A[5,1]
    Female = c(0,A[1:4,1]) + Female + c(0,A[8:11,1]) - MaleAGEP
    Male = c(0, A[1:4,1])
    Combined = cbind(Male,Female)

    SEmale = c(0, A[1:4,2])
    SEfemale = c(A[5,2], 
                sqrt(A[1:4,2]^2+A[8:11,2]^2+2*B[1:4]))
    SE = cbind(SEmale,SEfemale)
    Lower = exp(Combined - 2.576*SE)
    Upper = exp(Combined + 2.576*SE)
    Lower = apply(Lower, 2, percent)
    Upper = apply(Upper, 2, percent)
    Text = matrix(paste0("(", Lower, ", ", Upper, ")"), nrow = 5, ncol = 2, byrow = F)
    Text[1,1] = "Baseline"
    colnames(Text) = c("Male","Female")
    rownames(Text) = c("Married","Widowed","Divorced","Separated","Never Married")
    capt = "99% Confidence Intervals for Marital Status relative to Married Men for Hour model"
    myHeader = c(capt = 3)
    names(myHeader) = c(capt)
    A = knitr::kable(Text, "latex", booktabs = TRUE, linesep = "") %>%
        kable_styling(latex_options = c("HOLD_position")) %>%
        column_spec(1, bold = TRUE) %>%
        column_spec(1:3, border_left = TRUE) %>%
        add_header_above(header = myHeader, bold = TRUE)
    return(A)
}
    
PrintReplicationsTable <- function(Replications){
    UpperYear = Replications[[4]][1]+Replications[[6]]*2.576
    LowerYear = Replications[[4]][1]-Replications[[6]]*2.576
    UpperHour = Replications[[2]][1]+Replications[[5]]*2.576
    LowerHour = Replications[[2]][1]-Replications[[5]]*2.576

    Tab = matrix(c(Replications[[4]][1], Replications[[2]][1], LowerYear, LowerHour, UpperYear, UpperHour), nrow = 2, ncol = 3)
    colnames(Tab) = c("Estimate","99% Confidence Interval Lower Bound","Upper Bound")
    rownames(Tab) = c("Earnings per Year","Earnings per Hour")
    
    Tab = as.data.frame(Tab)
    capt = paste0("99% Confidence Interval for Adjusted R-squared")
    A = knitr::kable(Tab, "latex", booktabs = TRUE, linesep = "") %>%
      kable_styling(latex_options = c("HOLD_position")) %>%
      column_spec(2:4, border_left = TRUE)
    return(A)
}
requiredLibraries <- c("XML", "here", "tidyverse", "lubridate", "stats", "pracma", "IDPmisc", "ggplot2")
lapply(requiredLibraries, require, character.only = TRUE) # this is useful, as you can detect FALSES which means the library isn't loaded so you may as well abort

#which test to import and get the relevant files
testList    <- list.dirs(path = "DataFiles/", full.names = F, recursive = F)
folder      <- menu(testList, graphics = T, title = "Which test to analyse?") 
RawFiles    <- dir(path = str_c("DataFiles/",testList[folder]), pattern = "*.csv$", recursive = T, full.names = T)
TestTLA     <- c("Visit1", "Visit2")

#RawFiles   <- RawFiles[randi(length(RawFiles))]# change the number if you want a single file
#RawFiles   <- RawFiles[1]# change the number if you want a single file

sprintf("Found %i raw files in folder %s, and will collate them into %sdf.Rda", length(RawFiles), testList[folder], TestTLA[folder])

# Analysis outcomes based on https://iovs.arvojournals.org/article.aspx?articleid=2598502

# Group A or B - intervention or control but allocation unknown
# Session; 
#1 = Visit 1 pre cold pressor
#2 = Visit 1 during cold pressor
#3 = Visit 1 post cold pressor
#4 = Visit 2 pre cold pressor
#5 = Visit 2 during cold pressor
#6 = Visit 2 post cold pressor

#test


if (TestTLA[folder] == "Visit1") {  
  ImportFnc <- function(RawFiles){
    csvDF= read.csv(RawFiles)
    
    rawFilename     <- str_split_fixed(RawFiles[1], "/",3)
    ID <- rawFilename[3]
    ID <- as.numeric(str_extract(ID, "[0-9]+"))
    TestType     <- str_split_fixed(rawFilename, "_",3)
    TestType <- TestType[3,1]
    TestType <- as.character(str_extract(TestType, "[aA-zZ]+"))
 
    # Rename time coloumn
    csvDF <- csvDF %>%
      rename(Time = 4)
    
    tot <- nrow(csvDF)
    #Percentage inval data
    InVal <- nrow(filter(csvDF, LPMMV > 0, RPMMV > 0))
    percen_Inval <- (InVal / tot) *100
    #Percentage Left eye inval data
    InVal_L <- nrow(filter(csvDF, LPMMV == 1))
    percen_Inval_L <- (InVal_L / tot) *100
    #Percentage Right eye inval data
    InVal_R <- nrow(filter(csvDF, RPMMV > 0))
    percen_Inval_R <- (InVal_R / tot) *100
    # percentage pupil size over 10mm
    over_10 <- nrow(filter(csvDF, LPMM > 10, RPMM > 10))
    Pupil_10 <- (over_10 / tot) *100
    # percent pupil size under 1mm
    under_1 <- nrow(filter(csvDF, LPMM < 1, RPMM < 1))
    Pupil_1 <- (under_1 / tot) *100
    
    # cut down to what colomns needed
    #csvDF <- csvDF[c(1:3,7, 37:40)]
    
    
    #filter out pupil size values that are not "valid"
    csvDF <- filter(csvDF, LPMMV > 0, RPMMV > 0)

    #filter out pupil sizes that are below zero - should not exist but do
    csvDF <- filter(csvDF, LPMM > 1, RPMM > 1)
    
    #filter out pupil sizes that are above 10mm - should not exist but do (https://www.ncbi.nlm.nih.gov/books/NBK381/#:~:text=The%20normal%20pupil%20size%20in,opposite%20eye%20(consensual%20response))
    csvDF <- filter(csvDF, LPMM < 10, RPMM < 10)
    

    
    #mean of each pupil
    L_Mean <- mean(csvDF$LPMM, na.rm=T)
    R_Mean <- mean(csvDF$RPMM, na.rm=T)
    
    # MEan of the difference in pupil change over time
    L_Diff <- diff(csvDF$Time) / diff(csvDF$LPMM)
    L_Diff[is.infinite(L_Diff)] <- NA
    L_Diff_mean <- mean(na.omit(L_Diff))
    R_Diff <- diff(csvDF$Time) / diff(csvDF$RPMM)
    R_Diff[is.infinite(R_Diff)] <- NA
    R_Diff_mean <- mean(na.omit(R_Diff))
    

      
    data.frame(ID = ID,
               test_Type = TestType,
               both_Invalid = percen_Inval,
               L_Inval = percen_Inval_L,
               R_Inval = percen_Inval_R,
               under1mm = under_1, 
               over10mm = over_10,
               L_Mean_Pup_Diam = L_Mean,
               L_Mean_Pup_Diam_Change = L_Diff_mean,
               R_Mean_Pup_Diam = R_Mean,
               R_Mean_Pup_Diam_Change = R_Diff_mean
    )
    
  }
}else if (TestTLA[folder] == "Visit2") {  
  ImportFnc <- function(RawFiles){
    csvDF= read.csv(RawFiles)
    
    csvDF <- csvDF %>%
      rename(Time = 7)
    
    
    csvDF <- csvDF[c(1:3,7, 37:40)]
    
    #filter out pupil size values that are not "valid"
    csvDF <- filter(csvDF, LPMMV > 0, RPMMV > 0)
    
    #filter out pupil sizes that are below zero - should not exist but do
    csvDF <- filter(csvDF, LPMM > 0, RPMM > 0)
    
    #filter out pupil sizes that are above 10mm - should not exist but do (https://www.ncbi.nlm.nih.gov/books/NBK381/#:~:text=The%20normal%20pupil%20size%20in,opposite%20eye%20(consensual%20response))
    csvDF <- filter(csvDF, LPMM < 10, RPMM < 10)
    
    L_Mean <- mean(csvDF$LPMM)
    R_Mean <- mean(csvDF$RPMM)
    
    L_Diff <- diff(csvDF$Time) / diff(csvDF$LPMM)
    L_Diff_mean <- mean(NaRV.omit(L_Diff))
    
    R_Diff <- diff(csvDF$Time) / diff(csvDF$RPMM)
    R_Diff_mean <- mean(NaRV.omit(R_Diff))
    
    # Fast fourier transformation 
    #N <- nrow(csvDF)
    #csvDF$LPMM <- sin( 2*pi * csvDF$Time )
    #with( csvDF, plot( Time, LPMM, xlab = "t (sec)" ) )
    #
    #fdomain <- data.frame( f = seq( 0, N - 1 ) / ( N * diff( csvDF$Time[ 1:2 ] ) ) )
    #fdomain$mag <- Mod( fft( csvDF$LPMM) )
    #with( fdomain[ 1:11, ], plot( f, mag, xlab = "f (Hz)" ) ) # confirm that the 1sec per cycle sine shows up as a 1Hz spike
    #with( fdomain, plot( f, mag, xlab = "f (Hz)" ) ) # confirm symmetry
    #
    #magnitudes<- abs(fft_result)
    ## Find the frequency with the largest amplitude
    #max_index_sine <- which.max(magnitudes)
    #mag_Value <-  magnitudes[max_index_sine]
    
    
    
    data.frame(ID = first(csvDF$ID),
               GROUP= first(csvDF$GROUP), 
               SESSION = first(csvDF$SESSION),
               Time = max(csvDF$Time),
               L_Mean_Pup_Diam = L_Mean,
               L_Mean_Pup_Diam_Change = L_Diff_mean,
               R_Mean_Pup_Diam = R_Mean,
               R_Mean_Pup_Diam_Change = R_Diff_mean
    )
    
  }
}



collatedDF <- do.call(rbind.data.frame, lapply(RawFiles, ImportFnc))
save(collatedDF,file=sprintf("%sdf.Rda", TestTLA[folder]))
view(collatedDF)

write.csv(collatedDF, "collatedDF.csv")


require("kableExtra")
require("broom")

get_summary_stats(collatedDF, show = c("min", "max", "mean")) %>%
  kbl(caption = "Summary Data for Eye Tracking Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")

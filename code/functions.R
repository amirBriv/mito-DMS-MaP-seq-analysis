library(tidyr)
library(dplyr)
library(ggmsa)
library(circlize)
library(ggplot2)
library(zoo)
library(slider)
library(readr)
library(stringr)
library(pROC)
library(viridis)
library(ggrepel)
library(DescTools)
library(wesanderson)
library(mltools)
library(ROCR)
library(R4RNA)
library("scales")
#####Functions
read.popavg = function(path,GU=TRUE){
  x = read.table(path,sep = "\t",header = T)
  if(GU){
    return(x)
  }
  else{
    x_AC = get_AC(x)
    return(x)
  }
}

get_AC = function(df){
  sub_df = df %>%
    filter(Base =="A" | Base == "C") 
  return(sub_df)
}
set_GU_0 = function(df){
  df = df %>% 
    mutate(Mismatches...Deletions = case_when(Base == 'T' ~ 0,
                                              Base == 'G' ~0,
                                              Base == 'A' ~ Mismatches...Deletions,
                                              Base == 'C' ~ Mismatches...Deletions,
    ))
  
}

get_GU = function(df){
  sub_df = df %>%
    filter(Base =="G" | Base == "T") 
  return(sub_df)
}

# slideFunct <- function(data, window, step){
#   total <- length(data)
#   spots <- seq(from=1, to=(total-window), by=step)
#   result <- vector(length = length(spots))
#   pos = vector(length = length(spots))
#   for(i in 1:length(spots)){
#     result[i] <- mean(data[spots[i]:(spots[i]+window)])
#     pos[i] = i
#   }
#   return(result)
# }

make_message_table = function(pos_avg.list, WT_or_KO="WT", Message){
  if(WT_or_KO == "WT"){
    mWT1 = pos_avg.list$WT1 %>%
      filter(Gene == Message)
    mWT2 = pos_avg.list$WT2 %>%
      filter(Gene == Message)
    mWT1 = get_AC(mWT1)
    mWT2 = get_AC(mWT2)
    message_tibble = tibble(WT1 = mWT1$Mismatches...Deletions,
                            WT2 = mWT2$Mismatches...Deletions,
                            Base = mWT1$Base,
                            Position = mWT1$Position)
  }
  if(WT_or_KO == "KO"){
    mKO1 = pos_avg.list$KO1 %>%
      filter(Gene == Message)
    mKO2 = pos_avg.list$KO2 %>%
      filter(Gene == Message)
    mKO1 = get_AC(mKO1)
    mKO2 = get_AC(mKO2)
    message_tibble = tibble(KO1 = mKO1$Mismatches...Deletions,
                            KO2 = mKO2$Mismatches...Deletions,
                            Base = mKO1$Base,
                            Position = mKO1$Position)
  }
  return(message_tibble)
  
}


###Pop avg normailization from df
Normalize_Reactivities = function(df,threshold){
  
  ##Get ACs, only want to normalize to these reactivities
  df_AC = get_AC(df)
  #Put the reactivites in descending order (highest reactivities first)
  df_AC = df_AC[order(df_AC$Mismatches...Deletions, decreasing = T),]
  norm_thresh = as.integer(length(df_AC[,1])*threshold)
  #Get the median value of the top n bases where n is the threshold value provided (usually 10%) 
  median = median(df_AC$Mismatches...Deletions[1:norm_thresh])
  #Calculate Normalized Reactivites
  print(median)
  df$NormalizedReactivity = df$Mismatches...Deletions/median
  #Remove Ts and Gs
  df = df %>% 
    mutate(NormalizedReactivity = case_when(Base == 'T' ~ 0,
                                            Base == 'G' ~0,
                                            Base == 'A' ~ NormalizedReactivity,
                                            Base == 'C' ~ NormalizedReactivity,
    ))
  #Winsorization (fancy word for any of values greater than the median n% just set to 1)
  df = df %>% 
    mutate(NormalizedReactivity = case_when(NormalizedReactivity > 1 ~ 1,
                                            NormalizedReactivity == 1 ~ 1,
                                            NormalizedReactivity < 1 ~ NormalizedReactivity,
                                            NormalizedReactivity == NA ~0))
  #
  #df_norm = tibble(NormalizedReactivity = df$NormalizedReactivity
  return(df)
}

Calc_Gini = function(mu_list, na.rm=F) {
  if(na.rm){
    mu_list<-mu_list[!is.na(mu_list)]
  }
  n = length(mu_list)
  denom = 2 * n * sum(mu_list)
  numer = 0
  for(i in 1:n){
    for(j in 1:n){
      x_i = mu_list[i]
      x_j = mu_list[j]
      value = abs(x_i - x_j)
      numer = numer + value
    }
  }
  gini = numer / denom
  return(gini)
}

Read.CT = function(path){
  CT_RAW = read_tsv(path,show_col_types = FALSE)
  CT_df = tibble(index = c(),
                 base = c(),
                 paired = c())
  
  for(i in 1:length(CT_RAW[[1]])){
    CT_LINE = CT_RAW[i,]
    CT_LINE
    Line_Split = str_split(CT_LINE," ")
    Line_Split = Line_Split[[1]]
    Line_Split = Line_Split[nzchar(Line_Split)]
    Line_Split[1]
    index = Line_Split[1]
    base = Line_Split[2]
    paired = Line_Split[5]
    row_tibble = tibble(index = as.numeric(index),
                        Base = base,
                        paired = as.numeric(paired))
    CT_df = bind_rows(CT_df,row_tibble)
  }
  exit_row_index <- which(CT_df$Base == "ENERGY" & !is.na(CT_df$Base))
  if (length(exit_row_index) > 0) {
    CT_df = CT_df[1:(exit_row_index[1] - 1), ]
  }
  return(CT_df)
}



Build_ROC_df = function(CT_df,Pop_Avg_df,include_GU){
  ROC_df = tibble(Index = CT_df$index,
                  Base = Pop_Avg_df$Base,
                  Mu = Pop_Avg_df$Mismatches...Deletions,
                  Paired = CT_df$paired)
  ROC_df = ROC_df %>%
    mutate(Paired = case_when(Paired > 0 ~ 0,
                              Paired == 0 ~ 1))
  ROC_df
  if(!include_GU){
    ROC_df = subset(ROC_df, Base != "G")
    ROC_df = subset(ROC_df, Base != "T")
    ROC_df
  }
  else{
    return(ROC_df)
  }
}

Build_Paired_Table = function(CT_df){
  CT_df = CT_df %>%
    mutate(paired = case_when(paired > 0 ~ 0,
                              paired == 0 ~ 1))
}
compute_mFMI = function(CT_a, CT_b, include_GU){
  # if(!include_GU){
  #   CT_a = get_AC(CT_a)
  #   CT_b = get_AC(CT_b)
  # }
  P_common = 0
  P_unique_a = 0
  P_unique_b = 0
  u_unpaired = 0
  L = length(CT_a$paired)
  if(!(length(CT_a$index) == length(CT_b$index))){
    print("CT_Files need to be of same length!")
    break
  }
  
  for(i in 1:L){
    if(CT_a$paired[i] == 0 && CT_b$paired[i] == 0){
      u_unpaired = u_unpaired + 1
    }
    
    if(CT_a$paired[i]!=0 && CT_b$paired[i] != 0){
      if(CT_a$paired[i] != CT_b$paired[i]){
        P_unique_a = P_unique_a + 1
        P_unique_b = P_unique_b + 1
      }
      if(CT_a$paired[i] == CT_b$paired[i]){
        P_common = P_common + 1
      }
    }
    if(CT_a$paired[i] == 0 && CT_b$paired[i] != 0){
      P_unique_b = P_unique_b + 1
    }
    if(CT_a$paired[i] != 0 && CT_b$paired[i] == 0){
      P_unique_a = P_unique_a + 1
    }
  }
  u = u_unpaired/L
  FMI = P_common/(sqrt((P_common+P_unique_a)*(P_common+P_unique_b)))
  mFMI = u + ((1-u)*FMI)
  return(c(FMI,mFMI))
}

Gini_Slide = function(df, window_size, fill = TRUE, na.rm = F){
  df_AC = get_AC(df)
  Mus = df_AC$Mismatches...Deletions
  Positions = df_AC$Position
  Mu_Windows = slide(Mus, ~.x, .after  = window_size)
  Pos_Windows = slide(Positions, ~.x, .after  = window_size)
  
  Gini_Vect = c()
  Pos_Vect = c()
  if(na.rm){
    for(i in 1:length(Mu_Windows)){
      Gini_Vect[i] = Calc_Gini(Mu_Windows[[i]],na.rm = T)
      Pos_Vect[i] = Pos_Windows[[i]]
    }
  }
  else{
    for(i in 1:length(Mu_Windows)){
      Gini_Vect[i] = Calc_Gini(Mu_Windows[[i]],na.rm = F)
      Pos_Vect[i] = Pos_Windows[[i]]
    }
  }
  
  Window_df = tibble(Position = Pos_Vect,
                     Gini = Gini_Vect)
  Gini_Final_df = tibble(Position = df$Position)
  Gini_Final_df = left_join(Gini_Final_df,Window_df,"Position")
  if(fill){
    Gini_Fill_df = Gini_Final_df %>%
      fill(Gini, .direction = "down")
    return(Gini_Fill_df)
  }
  else{
    return(Gini_Final_df)
  }
}

Corr_Slide = function(df1,df2,window_size){
  if(!setequal(df1$Base,df2$Base) | length(df1$Position) != length(df2$Position)){
    print("Yo! These dfs have different Bases!")
  }
  else{
    df1_AC = get_AC(df1)
    df2_AC = get_AC(df2)
    tail(df1_AC)
    tail(df2_AC)
    df1_Mus = df1_AC$Mismatches...Deletions
    df2_Mus = df2_AC$Mismatches...Deletions
    
    df1_Positions = df1_AC$Position
    
    df1_Mu_Windows = slide(df1_Mus, ~.x, .after  = window_size)
    df2_Mu_Windows = slide(df2_Mus, ~.x, .after  = window_size)
    Pos_Windows = slide(df1_Positions, ~.x, .after  = window_size)
    Pos_Windows
    corr_vect = c()
    pos_vect = c()
    
    
    for(i in 1:length(df1_Mu_Windows)){
      corr_vect[i] = cor(df1_Mu_Windows[[i]], df2_Mu_Windows[[i]],method = "pearson", use = "pairwise.complete.obs")^2
      pos_vect[i] = median(Pos_Windows[[i]])
    }
    hold_tib = tibble(Position = pos_vect,
                      Corr = corr_vect )
    
    final_corr_tib = tibble(Position = df1$Position)
    final_corr_tib = left_join(final_corr_tib,hold_tib,"Position")
    final_corr_tib = final_corr_tib %>%
      fill(Corr, .direction = "down")
    return(final_corr_tib)
  }
}

Corr_Slide_R = function(df1,df2,window_size){
  if(!setequal(df1$Base,df2$Base) | length(df1$Position) != length(df2$Position)){
    print("Yo! These dfs have different Bases!")
  }
  else{
    df1_AC = get_AC(df1)
    df2_AC = get_AC(df2)
    tail(df1_AC)
    tail(df2_AC)
    df1_Mus = df1_AC$Mismatches...Deletions
    df2_Mus = df2_AC$Mismatches...Deletions
    
    df1_Positions = df1_AC$Position
    
    df1_Mu_Windows = slide(df1_Mus, ~.x, .after  = window_size)
    df2_Mu_Windows = slide(df2_Mus, ~.x, .after  = window_size)
    Pos_Windows = slide(df1_Positions, ~.x, .after  = window_size)
    Pos_Windows
    corr_vect = c()
    pos_vect = c()
    
    
    for(i in 1:length(df1_Mu_Windows)){
      corr_vect[i] = cor(df1_Mu_Windows[[i]], df2_Mu_Windows[[i]],method = "spearman", use = "pairwise.complete.obs")
      pos_vect[i] = median(Pos_Windows[[i]])
    }
    hold_tib = tibble(Position = pos_vect,
                      Corr = corr_vect )
    
    final_corr_tib = tibble(Position = df1$Position)
    final_corr_tib = left_join(final_corr_tib,hold_tib,"Position")
    final_corr_tib = final_corr_tib %>%
      fill(Corr, .direction = "down")
    return(final_corr_tib)
  }
}


AUROC_Slide = function(mu_df, CT_File, window_size) {
  roc_df = Build_ROC_df(CT_df = CT_File, Pop_Avg_df = mu_df, include_GU = T)
  start_pos = 1
  # Adjusting the size of the vectors to match the possible window positions
  auc_vect = numeric(length(roc_df$Mu) - window_size + 1)
  pos_vect = numeric(length(roc_df$Mu) - window_size + 1)
  
  for (i in 1:(length(roc_df$Mu) - window_size + 1))  {
    temp_ROC = roc_df[start_pos:(start_pos + window_size - 1), ]
    temp_ROC_AC = get_AC(temp_ROC)
    temp_ROC_AC = na.omit(temp_ROC_AC)
    auc_vect[i] = auc_roc(temp_ROC_AC$Mu, temp_ROC_AC$Paired)
    pos_vect[i] = start_pos + (window_size/2) - 1
    start_pos = start_pos + 1
  }
  temp_tib = tibble(Position = pos_vect, AUC = auc_vect)
  return(temp_tib)
}


Make_Arc = function(path){
  R4RNA_Obj = readConnect(path)
  R4RNA_Helix = expandHelix(R4RNA_Obj)
  R4RNA_Helix$col = "blue"
  plotHelix(R4RNA_Helix,line = T)
}

Make_Double_Arc = function(path1,path2){
  struct_WT = R4RNA::readConnect(path1)
  struct_Helix_WT = expandHelix(struct_WT)
  struct_Helix_WT$col = "purple"
  struct_KO = readConnect(file = path2)
  struct_Helix_KO = expandHelix(struct_KO)
  struct_Helix_KO$col = "lightblue"
  plotDoubleHelix(struct_Helix_WT,struct_Helix_KO,line = T)
}


set_outliers_na <- function(data_list, threshold) {
  #Make new list to return
  # Loop through each data frame in the list
  na_list = list()
  for (i in 1:length(data_list)) {
    # Get the current data frame
    df <- data_list[[i]]
    na_list[[i]] = df[df$Mismatches...Deletions > threshold,]
  }
  #get all positions (unique) with outliers
  na_pos = unique(bind_rows(na_list)$Position)
  # Set the rows corresponding to the outliers to NA
  for (i in 1:length(data_list)) {
    # Get the current data frame
    df <- data_list[[i]]
    # Set the rows corresponding to the outliers to NA
    df[na_pos,"Mismatches...Deletions"] = NA
    # Replace the original data frame in the list with the modified data frame
    data_list[[i]] <- df
  }
  return(data_list)
}

Convert_Neg_Coords = function(df){
  df=df[order(nrow(df):1),]
  df$Position = 1:length(df$Position)
  return(df)
}

write.normalized_reacts = function(path,name){
  df1 = read.popavg(path)
  df1_Norm = Normalize_Reactivities(df1,0.1)
  df1_tibble = tibble(Position = df1_Norm$Position,
                      React = df1_Norm$NormalizedReactivity)
  struct_path = sub('/[^/]*$', '', path)
  file_name = paste(struct_path,"/",name,".csv", sep= "")
  write.csv(df1_tibble,file_name,quote = F)
}


#############CLUSTER FUNCTIONS###############

Get_Whole_Cluster_Correlation = function(K2_df) {
  good_ratio_df <- K2_df %>%
    filter(ratio == 1) %>%
    mutate(chunk = (row_number() - 1) %/% 100)  # Define chunks
  good_ratio_df_AC = get_AC(good_ratio_df)
  # Compute correlation for each chunk
  correlation_results <- good_ratio_df_AC %>%
    group_by(chunk) %>%
    summarize(correlation = cor(Cluster_1, Cluster_2, method = 'pearson', use = "complete.obs")^2, .groups = 'drop')
  good_ratio_df_AC = left_join(good_ratio_df_AC,correlation_results)
  good_ratio_df = left_join(good_ratio_df,good_ratio_df_AC)
  good_ratio_df = good_ratio_df %>%
    fill(correlation, .direction = "down")
  K2_df = left_join(K2_df,good_ratio_df)

  return(K2_df)
 }

####Circos Functions###########
Convert_Circos_Coords = function(df){
  R_Loop = df %>%
    filter(Position < 577)
  Coding_Region = df %>%
    filter(Position > 576)%>%
    filter(Position < 16569)
  tail(R_Loop)
  head(Coding_Region)
  Circos_df = bind_rows(Coding_Region,R_Loop)
  Circos_df$Position = 1:length(Circos_df$Position)
  return(Circos_df)
}



Init_Circos_Mito = function(){
  pal <- wes_palette("Darjeeling1", 5, type = "discrete")
  circos.csv = read.csv("../data/Circos_Table.csv")
  mito_length = "16568"
  ticks = seq(0,17001,1000)
  df = data.frame(start = 0, end = 16568)
  rownames(df) = c("Mitochondria")

  #pal = pal_aaas()(10)
  circos.par("start.degree" = 90,"canvas.xlim" = c(-1.5,1.5))
  circos.par("gap.degree" = 0)
  circos.initialize(xlim = df)
  circos.par("track.height" = 0.2)
  circos.track(ylim = c(1,10))
  circos.axis(h = "top",
              major.at = ticks,
              labels.facing = "clockwise",
              major.tick.length = 2,
              labels.cex = 1)

  for(i in 1:length(circos.csv$Gene)){
    if(circos.csv$Strand[i] == '+'){
      if(circos.csv$Type[i] == "tRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5,col = pal[1])

      }
      if(circos.csv$Type[i] == "mRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5, col = pal[5])

      }
      if(circos.csv$Type[i] == "rRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5,col = "grey")

      }
      if(circos.csv$Type[i] == "NC"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 1,col = "grey")
      }
    }
    if(circos.csv$Strand[i] == '-'){
      if(circos.csv$Type[i] == 'tRNA'){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 5,ybottom = 1,col = pal[1])

      }

      if(circos.csv$Type[i] == 'mRNA'){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 5,ybottom = 1,col = pal[5])
      }
    }
    if(circos.csv$End[i] < 8630){
      circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 25, circos.csv$Gene[i], facing = "clockwise", cex = 0.6, niceFacing = T)
    }else{
      circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 25, circos.csv$Gene[i], facing = "reverse.clockwise", cex = 0.6)

    }
  }

}

Convert_Neg_Coords = function(df){
  df=df[order(nrow(df):1),]
  df$Position = 1:length(df$Position)
  return(df)
}

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

DMS_React_Mean_Sliding = function(df1,window_size,AC=T,GU=F){
  if(AC && GU || !AC && !GU){
    print("Calculating Avg DMS reactivity on ACs and GUs")
    df1_Signal = df1
    Positions = df1$Position
  }
  if(AC){
    df1_Signal = get_AC(df1)
    Positions = df1_Signal$Position
  }
  if(GU){
    df1_Signal = get_GU(df1)
    Positions = df1_Signal$Position
  }
  Mus = df1_Signal$Mismatches...Deletions
  Mu_Windows = slide(Mus, ~.x, .after  = window_size)
  Pos_Windows = slide(Positions, ~.x, .after  = window_size)

  Mu_Mean_Vect = c()
  Position_Vector = c()
  length(Mu_Windows)
  for(i in 1:length(Mu_Windows)){

    Mu_Mean_Vect[i] = mean(Mu_Windows[[i]],na.rm=TRUE)
    Position_Vector[i] = median(Pos_Windows[[i]])
  }
  df_windows = tibble(Position = Position_Vector,
                      DMS_Mean = Mu_Mean_Vect)
  df_final = tibble(Position = df1$Position)
  df_final = left_join(df_final,df_windows)
  df_final = df_final %>%
    fill(DMS_Mean, .direction = "down")
  df_final
  return(df_final)
}






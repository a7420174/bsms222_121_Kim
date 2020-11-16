library(tidyverse)
library(readxl)
S1_1 <- read_excel("aav8130_Data_S1.xlsx")
S4_4 <- read_excel("aav8130_Data_S4.xls", sheet = 4)
S1_1 <- S1_1[-1,]
colnames(S1_1) <- c(colnames(S1_1)[1:7],'A','Bv','Bnv','C','D',colnames(S1_1)[13])
S1_1 <- S1_1 %>% filter(!(A=='NA')) %>% mutate('Patient ID' = as.character(`Patient ID`), A = as.numeric(A))


ind_DEGs <- as.data.frame(sapply(S1_1$`Patient ID`, 
                     function(ID) str_count(S4_4$`ASD patients with DE signal`, ID)))
S4_4 <- cbind(S4_4, ind_DEGs)

DEGs_num <- S4_4 %>% group_by(`Cell type`) %>% summarise_at(vars(`5278`:`4849`), sum)
tDEGs_num <- t(DEGs_num[,-(1)])
rownames(tDEGs_num) <- NULL
colnames(tDEGs_num) <- DEGs_num$`Cell type`
df <- data.frame(PatientID = S1_1$`Patient ID`, A = S1_1$A, tDEGs_num)
df %>% ggplot(aes(A, , col = PatientID)) + geom_point()
cor(df$A, df[,-(1:2)], method = "spearman")


df1 <- df %>% filter(PatientID != "5978")
cor(df1$A, df1[,-(1:2)], method = "pearson")
df1 %>% ggplot(aes(A, IN.VIP, col = PatientID)) + geom_point()

DEGs_num2 <- S4_4 %>% group_by(`Cell type`, `Direction of change (ASD/Control)`) %>% summarise_at(vars(`5278`:`4849`), sum)
tDEGs_num2 <- t(DEGs_num2[,-(1:2)])
rownames(tDEGs_num2) <- NULL
colnames(tDEGs_num2) <- paste(DEGs_num2$`Cell type`, DEGs_num2$`Direction of change (ASD/Control)`, sep = "_")
df2 <- data.frame(PatientID = S1_1$`Patient ID`, A = S1_1$A, tDEGs_num2)
df3 <- df2 %>% filter(PatientID != "5978")
cor(df3$A, df3[,-(1:2)], method = "pearson")
df3 %>% ggplot(aes(A, L4_UP, col = PatientID)) + geom_point()

S1_1 <- S1_1 %>% filter(!(A=='NA'&B1=='NA'&B2=='NA'&C=='NA'&D=='NA')) %>%
  mutate(rankA = rank(as.numeric(A), na.last = "keep"), 
         rankB1 = rank(as.numeric(Bv), na.last = "keep"), 
         rankB2 = rank(as.numeric(Bv), na.last = "keep"), 
         rankC = rank(as.numeric(C), na.last = "keep"), 
         rankD = rank(as.numeric(D), na.last = "keep"))
rankB <- vector()
  for (i in 1:nrow(S1)){
    if (is.na(S1$rankB1[i])){
      rankB[i] = S1$rankB2[i]
    }
    else if (is.na(S1$rankB2[i])){
      rankB[i] = S1$rankB1[i]
    }
    else {
      rankB[i] = (S1$rankB1[i]+S1$rankB2[i])/2
    }
  }
S1 <- S1 %>% mutate(rankB = rankB,
                    combinedscore = rankA + rankB + rankC + rankD)

y <- function(index) {
  ids <- str_extract_all(S4_4$'ASD patients with DE signal', '[0-9]+')
  ids[[index]] %in% c("5278", "5278", "5419", "5565", "5864", "5939", "5945", "6033")
}
sapply(1:nrow(S4_4), y)
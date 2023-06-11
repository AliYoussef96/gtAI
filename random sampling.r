library(infer)
library(dplyr)

library(openxlsx)
library(stats)
library(ggplot2)



df <- read.xlsx("Table S4 (Eukarya, fungal).xlsx", sheet = "Aspergillus fumigatus")

sample.size <- 0.25 * nrow(df)



df.r <- df %>% rep_sample_n(size = sample.size, reps = 1000 , replace = T) %>% 
  summarise(rgc = cor.test(gtAI, CAI, method = "spearman")$estimate,  
            pgc = cor.test(gtAI, CAI, method = "spearman")$p.value, 
            rsc = cor.test(StAI, CAI, method = "spearman")$estimate, 
            psc = cor.test(StAI, CAI, method = "spearman")$p.value,

            roc = cor.test(tAI, CAI, method = "spearman")$estimate,
            poc = cor.test(tAI, CAI, method = "spearman")$p.value)

df.r$gcfdr <- p.adjust(df.r$pgc , method = "fdr")
df.r$scfdr <- p.adjust(df.r$psc , method = "fdr")
df.r$ocfdr <- p.adjust(df.r$poc , method = "fdr")


df.r <- df.r[df.r$gcfdr < 0.05 && df.r$scfdr < 0.05 && df.r$poc < 0.05,]

df.r.plot <- df.r[,c(1,2,4,6)]
df.r.plot <- reshape2::melt(df.r.plot,id.vars = "replicate" )

p <- ggplot(data=df.r.plot, aes(x= replicate, y=value, color=variable)) + theme_classic()+
  geom_line(size = 1, alpha = 0.8) + scale_color_manual(labels = c("gtAI Vs. CAI", "stAI Vs. CAI", "otAI Vs. CAI"), values = c("#C0392B", "#7F8C8D", "#17A589")) + 
  theme(legend.title = element_blank()) + theme(text=element_text(size=10, face = "bold", colour = "black"),
                                                 axis.text = element_text(colour = "black", size =10 )) +
  ylab("rho values") + xlab("Replicates") + theme(legend.position="top") 

p 

#=====================================================================
# Yingying Dong.2020-8-25.Calculate weught.
#=====================================================================
setwd('/home/hp/weight/')
spe_name <- file("all_name.txt", open = "r")
spe=readLines(spe_name,n=1)
rscu_path = '/home/hp/RSCU/'
while( length(spe) != 0 ) {
  print(spe)
  write.table(spe,file = 'spe_log.txt',sep = '\n',quote = F,row.names = F,col.names = F)
  df = read.table(paste0(rscu_path,spe),header = T,sep = '\t',quote = "")
  df$Weights <- ave(df$RSCU,df$AA,FUN=function(x) x/max(x)) # calculate weight
  df = df[,-c(1,3,4,5)]
  df = df[-c(30,61,62,63,64),]
  write.table(df,file = paste0(spe,'_weight.txt'),sep = '\t',quote = F,row.names = F,col.names = F) 
  spe=readLines(spe_name,n=1)
}
close(spe_name)




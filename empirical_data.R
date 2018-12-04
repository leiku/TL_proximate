
############# get clean data ###########

dat.aphid<-readRDS("Data/Data_aphids.RDS")
dat.northsea<-readRDS("Data/Data_northsea.RDS")
dat.calcofi<-readRDS("Data/Data_chla.RDS")
name.data<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")

D<-vector("list",6)
names(D)<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
D[[1]]<-dat.aphid
D[[2]]<-dat.northsea
for(i in 3:6){D[[i]]<-dat.calcofi}

# format calcofi data in different groups
site<-as.numeric(row.names(dat.calcofi[[1]]))
group<-cut(site %% 1000,c(0,55,70,90,120),c(1,2,3,4))
for(gr in 1:4){
  ii<-which(group==gr)
  for(d in 1:10){
    D[[2+gr]][[d]]<-dat.calcofi[[d]][ii,]
  }
}
saveRDS(D,"data_original.RDS")

# species names
names1<-data.frame(dataset=c(rep("aphid",20),rep("northsea",22),rep("calcofi",10)), index=c(1:20,1:22,1:10),
                   names=c(names(D[[1]]),names(D[[2]]),names(D[[3]])))
write.csv(names1,"species_names.csv",row.names = F)

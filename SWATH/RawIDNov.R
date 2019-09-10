RID=read.table("pppb1566.txt", as.is = T, sep="\t") #original file
colnames(RID)[c(1,2)]=c("Raw_ID", "PPPA_ID")
RID$Raw_ID=toupper(RID$Raw_ID)
RID$PPPA_ID=toupper(RID$PPPA_ID)
RID$CPP_ID="NA" #CPP_ID
RID$PPP_ID="NA" #PPP_ID
RID$PAIR_ID="NA"
RID$MS_ID="NA" #MS_ID
RID$USZ_ID="NA" #USZ_ID
RID$Tissue="NA" # Tissue, C (ctrl), A (AQUA), N (normal), T (Cancer)
RID$Batch="NA" #Batch id
RID$Date="NA" #date
RID$Year="NA"
RID$Month="NA"
RID$Day="NA"
###extract CPP ID from raw id
pc_idx=c()
for(i in 1:nrow(RID))
{
  RID$Raw_ID[i]=sub("__", "_", RID$Raw_ID[i])
  RID$CPP_ID[i]=unlist(strsplit(RID$Raw_ID[i], split = "_"))[4]
  
  strL=unlist(strsplit(RID$Raw_ID[i], split = "_"))
  if(substr(strL[2], 1, 2) == "PC")
  {
    RID$MS_ID[i]=unlist(strsplit(RID$Raw_ID[i], split = "_"))[2]
    RID$Date[i]=unlist(strsplit(RID$Raw_ID[i], split = "_"))[3]
  } else {
    idxPC=which(substr(RID$MS_ID, 1, 2) != "PC")
    pc_idx=c(pc_idx, i)
    RID$MS_ID[i]=unlist(strsplit(RID$Raw_ID[i], split = "_"))[3]
    RID$Date[i]=unlist(strsplit(RID$Raw_ID[i], split = "_"))[2]
  }
}

for(i in 1:nrow(RID))
{
  RID$Year[i] = substr(RID$Date[i], 1, 2)
  RID$Month[i] = substr(RID$Date[i], 3, 4)
  RID$Day[i] = substr(RID$Date[i], 5, 6)
}

#idx

####trouble shooting 0: find disorded PC and date
write.table(RID$Raw_ID[pc_idx], "pc_error.txt", row.names = F, col.names = F, quote = F)
####trouble shooting 1: find incorrect sw id
sw_idx=grep("SW$", RID$Raw_ID, perl = T, invert=T)
write.table(RID$Raw_ID[sw_idx], "trouble_shooting_sw.txt", row.names = F, col.names = F, quote = F)

####trouble shooting 2: fix incorrect cp id
cp_idx=grep("^CP\\d+", RID$CPP_ID, perl=T)
write.table(RID[cp_idx,1], "trouble_shooting_cp.txt", row.names = F, col.names = F, quote=F)
RID$CPP_ID[cp_idx]=sub("CP", "CPP", RID$CPP_ID[cp_idx])
RID$Raw_ID[cp_idx]=sub("CP", "CPP", RID$Raw_ID[cp_idx])


##find batch id
for (i in 1:34)
{
  if (file.exists(paste0("./batch/", "batch", i, ".txt")))
  {
    bt=read.table(paste0("./batch/", "batch", i, ".txt"), as.is=T)
    for (j in 1:nrow(bt))
    {
      s1=unlist(strsplit(bt[j,], split="\\s+", perl = T))
      idx=(which(RID$CPP_ID == s1[1]))
      if (length(idx) != 0)
      {
        RID$Batch[idx] = i
        RID$PPP_ID[idx] = toupper(s1[2])
      }
    }
  }
}


##########
vID=read.csv("VarID.csv", header = T)
vID$cgb_pid0 = "NA"
vID$cgb_tissue = "NA"
for(i in 1:nrow(vID))
{
  id=unlist(strsplit(as.character(vID$rID[i]), split=" "))
  if(length(id)>1)
  {
    vID$cgb_pid0[i] = paste0(id[1], "_", id[2])
    vID$cgb_tissue[i]=id[3]
  }
}

##find tissue tag, usz_id
for(i in 1:nrow(RID))
{
  idx=which(as.character(vID$number.of.progressed.samples) == RID$CPP_ID[i])
  if(length(idx) == 1)
  {
    RID$USZ_ID[i] = vID$cgb_pid0[idx]
    RID$Tissue[i] = vID$cgb_tissue[idx]
  }
}


##find batch id
batchCNT = 0
for (i in 1:34)
{
  if (file.exists(paste0("./batch/", "batch", i, ".txt")))
  {
    batchCNT = batchCNT + 1
    bt=read.table(paste0("./batch/", "batch", i, ".txt"), as.is=T)
    for (j in 1:nrow(bt))
    {
      s1=unlist(strsplit(bt[j,], split="\\s+", perl = T))
      idx=(which(RID$CPP_ID == s1[1]))
      if (length(idx) != 0)
      {
        RID$Batch[idx] = batchCNT
        RID$PPP_ID[idx] = toupper(s1[2])
      }
    }
  }
}

#tissue for ctrl, aqua
idx_c=grep("ctrl", RID$PPP_ID, ignore.case = T)
RID$Tissue[idx_c] = "C"
#there are typo, so using "crtl" to match
idx_c1=grep("crtl", RID$PPP_ID, ignore.case = T)
RID$Tissue[idx_c1] = "C"

idx_a=grep("AQUA", RID$PPP_ID, ignore.case = T)
RID$Tissue[idx_a] = "A"
idx_a2=grep("POOL", RID$CPP_ID, ignore.case = T)
RID$Tissue[idx_a2] = "A"
###2UG samples are aqua
idx_2UG=grep("2UG", RID$Raw_ID, perl = T)
RID$Tissue[idx_2UG] = "A"
RID$Batch[idx_2UG] = 0


idx_CTRL0=grep("GUOT_PC1_170206_PCPOOL18_SW", RID$Raw_ID, perl = T)
RID$Batch[idx_CTRL0] = 18
idx_CTRL1=grep("GUOT_PC2_170220_PCPOOL18_SW", RID$Raw_ID, perl = T)
RID$Batch[idx_CTRL1] = 18
idx_CTRL2=grep("GUOT_PC3_170213_PCPOOL20B_SW", RID$Raw_ID, perl = T)
RID$Batch[idx_CTRL2] = 20
idx_CTRL3=grep("GUOT_PC4_170301_PCPOOL20_SW", RID$Raw_ID, perl = T)
RID$Batch[idx_CTRL3] = 20

##
repID=grep("\\D$", RID$CPP_ID, perl=T)
#B
repIDv=c()
for(i in 1:length(repID))
{
  Cid=substr(RID$CPP_ID[repID[i]], nchar(RID$CPP_ID[repID[i]]), nchar(RID$CPP_ID[repID[i]]))
  if(Cid == "B" || Cid == "C")
  {
    Bid=substr(RID$CPP_ID[repID[i]], 1, nchar(RID$CPP_ID[repID[i]])-1)
    idx=which(RID$CPP_ID == Bid)
    if (length(idx) == 1)
    {
      RID$USZ_ID[repID[i]] = RID$USZ_ID[idx]
      RID$Tissue[repID[i]] = RID$Tissue[idx]
      RID$Batch[repID[i]] = RID$Batch[idx]
    } else if(length(idx) > 1) {
      RID$USZ_ID[repID[i]] = RID$USZ_ID[idx[1]]
      RID$Tissue[repID[i]] = RID$Tissue[idx[1]]
      RID$Batch[repID[i]] = RID$Batch[idx[1]]
      repIDv=c(repIDv, i)
    }
  }
}

# #find batch friend by matching date and PC#ID
# BatNA=which(RID$Batch == "NA")
# RID_BatNA=RID[BatNA,]
# RID_Bat=RID[-BatNA,]
# bID=array(0, dim=length(BatNA))
# for(i in 1:length(BatNA))
# {
#   idx=which(RID_Bat$MS_ID == RID_BatNA$MS_ID[i] & RID_Bat$Date == RID_BatNA$Date[i])
#   bID[i]=length(table(RID_Bat$Batch[idx]))
#   if(bID[i] == 1)
#   {
#     RID$Batch[BatNA[i]] = RID_Bat$Batch[idx[1]]
#   }
# }

USZNA=which(RID$USZ_ID == "NA")
TisNA=which(RID$Tissue == "NA")
BatNA1=which(RID$Batch == "NA")
RID$Batch[BatNA1] = "TBD"
write.table(RID$Raw_ID[BatNA1], "BatchNA1.txt", row.names = F, col.names = F, quote = F)


PPPA_F=read.table("PPP.txt", as.is = T, header = T)
for(i in 1:nrow(PPPA_F))
{
  idx=which(RID$PPPA_ID == PPPA_F$PPPA_ID[i])
  if(length(idx) == 1)
  {
    RID$PAIR_ID[idx] = PPPA_F$PPP_ID[i]
  }
  else if(length(idx) > 1)
  {
    print(idx)
  }
}
RID$PPP_ID=RID$PAIR_ID
RID=RID[, -which(colnames(RID) == "PAIR_ID")]
#result
write.csv(RID, "PPPB1_lineup.csv", quote = F, row.names = F)

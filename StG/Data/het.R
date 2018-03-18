raw=read.table("het1.raw", as.is = T, header = T)
freq=colMeans(raw[,7:ncol(raw)])/2
EHet=sum(2*freq*(1-freq))
OHet=array(0, nrow(raw))
for(i in 1:nrow(raw))
{
  OHet[i] = length(which(raw[i, 7:ncol(raw)] == 1))
}
plot(OHet)

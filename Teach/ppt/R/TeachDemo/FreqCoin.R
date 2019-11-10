rp=1000
p=100
obs=array(0, rp)
for(i in 1:rp) {
  obs[i]=mean(rbinom(p, 1, 0.5))
}
hist(main="Coins", obs, freq=F, breaks = 20)

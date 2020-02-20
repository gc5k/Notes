#### input: bv, disease 0/1, sex
#### output: matrix(1:10, c(decile, or.lo, or, or.hi))
getDecileOrs = function(bv, phen, sex=NA) {
  goodphen = which(phen %in% 0:1)
  bv = bv[goodphen]
  phen = phen[goodphen]
  sex = sex[goodphen]
  mat = matrix(NA, 4, 10 )
  d = bv
  q = quantile(d, seq(0, 1, len=11))
  dec = list()
  for(i in 1:(length(q)-1)) { dec[[i]] = which(d >= q[i] & d < q[i+1])}
  for(decile in 1:10) {
    d.top = dec[[decile]]
    d.bot = dec[[1]]
    tb = c(d.top, d.bot)
    if(!all(is.na(sex)))
      df = data.frame(phen=phen[tb], group=c(rep(1, length(d.top)), rep(0, length(d.bot))), sex=sex[tb])
    else
      df = data.frame(phen=phen[tb], group=c(rep(1, length(d.top)), rep(0, length(d.bot))))
    logm = glm(phen ~ ., data=df, family=binomial(logit))
    #logm = glm(phen ~ group, data=df, family=binomial(logit))
    est = summary(logm)$coefficients[2,1]
    std.err = as.numeric(summary(logm)$coefficients[2,2])
    
    #or.raw = (sum(fam[d.top, 6] == 2) / sum(fam[d.top, 6] == 1)) / (sum(fam[d.bot, 6] == 2) / sum(fam[d.bot, 6] == 1))
    or.adj = as.numeric(exp(est))
    or.adj.lo = exp(est - 1.96 * std.err)
    or.adj.hi = exp(est + 1.96 * std.err)
    pred=1
    mat[1, decile + (pred-1)*10] = decile
    mat[2, decile + (pred-1)*10] = or.adj.lo
    mat[3, decile + (pred-1)*10] = or.adj
    mat[4, decile + (pred-1)*10] = or.adj.hi
    #mat[5, decile + (pred-1)*10] = or.raw
  }
  out = data.frame(t(mat))
  names(out) = c('decile', 'OR.low', 'OR', 'OR.high')
  out
}



plot.deciles.2 = function(df, name='profile score decile ORs', size=12) {
  
  require(ggplot2)
  require(grid)
  
  th2 = theme(panel.background = element_rect(fill = "white", colour='black'), text=element_text(size=size), legend.title=element_blank(), axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-2.4), plot.margin=unit(c(5,5,15,12),"mm"), panel.grid=element_blank(), legend.position="top")
  limits = aes(ymax=OR.high, ymin=OR.low)
  ar=.5
  ymax=10
  ro = 3
  p = ggplot(df, aes(x=decile, y=OR, colour=class )) + geom_pointrange(limits, size=1) + xlab('decile compared to lowest') + ylab('odds ratio (prediction accuracy)') + scale_x_continuous(breaks=1:10) #+ th2 + 
    #theme(aspect.ratio=ar) + theme(legend.key.size=unit(8, 'mm'), legend.position=c(.7,.3), legend.key.height=unit(3,"line")) +
    #geom_hline(yintercept=1) + coord_cartesian(ylim=c(0,ymax))
  p
  
}



n=1e3
data = getDecileOrs(rnorm(n), sample(0:1, n, replace=T))
plot.deciles.2(data, size=24)








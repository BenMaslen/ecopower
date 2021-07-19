library(mvabund)
data(spider)
spiddat = mvabund(spider$abund)

X = data.frame(spider$x)
X$Treatment2 = rep(c("A","B"),each=14)
X$Treatment4 = rep(c("A","B","C","D"),each=7)
X$month = rep(c("Jan","Feb","Mar","Apr"),each=7)

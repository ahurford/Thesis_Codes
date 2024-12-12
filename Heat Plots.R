require(deSolve)
require(ggplot2)
require(patchwork)
require(ggpubr)
cols <- c("#D55E00","#E69F00","#F0E442","#66CC99","#009E73","white", "#CC79A7","#999999","#56B4E9")

##### Parameters
# constraint on community member isolation
C1max <- 1000
# constraint on traveller isolation
C2max <- 30
# transmission rate
beta<-.0002
# importation rate
theta<-1
mu<-0.334
# relative transmissibility of travellers
c<-1
T<-1000

# This is designed to trigger termination of the simulation when I1>0.5
rootfun <- function(t, y, parameters){
  I1 = y[2]
  # However, with importations, it doesn't make sense to have this as
  # the endpoint of the simulation, so instead I just set Imin = 0 and have
  # the simulation end at T=100.
  y1 = 1 - I1
  return(y1)
}

SI<-function(t, y, parameters){
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  C1 = y[4]
  C2 = y[5]

  if(C1 > C1max){
    u1 = 0
  }
  if(C2 > C2max){
    u2 = 0
  }

  dS = - beta*S*(I1 + c*I2)
  dI1 = beta*S*(I1 + c*I2) - mu*I1 - u1*I1
  dI2 = theta - 2*mu*I2 - u2*I2
  dC1 = u1*I1
  dC2 = u2*I2
  return(list(c(dS,dI1,dI2,dC1,dC2)))
}
yinc = .025
xinc = .025
stratMx = NULL
u1vec <- seq(2,0,-xinc)
u2vec <- seq(3,0,-yinc)

# Pre-allocation
max_elim = rep(max(u2vec)+yinc,length(u1vec))
max_supp = max_elim
max_supp_circ = max_elim
max_circ_mit = rep(min(u2vec),length(u1vec))

for(i in seq(1,length(u1vec))){
  u1 = u1vec[i]
  if(i%%10==0){
    print(u1)
  }
  for(j in seq(1,length(u2vec))){
    u2 = u2vec[j]
  out <- ode(y = c(S=5000,I1=10,I2=1,C1=0,C2=0), parms = NULL, times = seq(0,T,1), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun, method = "radau")
  out<- data.frame(out)
 if(max(out$I1)>10 & max(out$C1)<=C1max & max(out$C2)>C2max & max(out$time)<T){
  outcome = "[mit, circ]"
  max_supp_circ[i] = u2+yinc/2 
} else if(max(out$I1)>10 & max(out$C1)>C1max & max(out$C2)<=C2max& max(out$time)<T){
  outcome = "[circ, cont]"
} else if(max(out$I1)>10 & max(out$C1)<C1max & max(out$C2)<C2max& max(out$time)<T){
  outcome = "[mit, cont]"
  max_supp[i] = u2+yinc/2
} else if(max(out$time)>=T){
  outcome = "no end"
  print("no end")
} else if(max(out$I1)>10 & max(out$C1)>C1max & max(out$C2)>C2max& max(out$time)<T){
  outcome = "[circ, circ]"
  max_circ_mit[i]=u2+yinc/2
}
  else if(max(out$I1)<=10 & max(out$time)<T & max(out$C1)<=C1max & max(out$C2)<=C2max){
    outcome = "[elim, elim]"
    max_elim[i] = u2+yinc/2
  }


J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
stratMx = rbind(stratMx, data.frame(u1 = u1, u2 = u2, strategy = outcome, J=J, T = max(out$time)))
}}

elim = data.frame(u1vec, max_elim, max_supp, max_supp_circ)
g1=ggplot() +
  geom_raster(data = stratMx,aes(x=u1, y=u2, fill= strategy),hjust = 1, vjust=1)+
  geom_line(data=elim,aes(x=u1vec+1.5*xinc,y=max_elim-2*yinc),col="white")+
  geom_line(data=elim,aes(x=u1vec,y=max_supp),col="white")+
  geom_line(data=elim,aes(x=u1vec,y=max_supp_circ),col="white")+
  geom_line(data=elim,aes(x=u1vec,y=max_circ_mit),col="white")+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  ggtitle("Low importations")+
  scale_fill_manual(values=cols, breaks = c("[elim, elim]", "[mit, cont]", "[mit, circ]", "[circ, cont]", "[circ, circ]"))+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size

theta <- 2
stratMx2 = NULL
# Pre-allocation
max_elim2 = rep(max(u2vec)+yinc,length(u1vec))
max_supp2 = max_elim2
max_supp_circ2 = max_elim2
max_circ_mit2 = max_circ_mit

for(i in seq(1,length(u1vec))){
  u1 = u1vec[i]
  if(i%%10==0){
    print(u1)
  }
  for(j in seq(1,length(u2vec))){
    u2 = u2vec[j]
    out <- ode(y = c(S=5000,I1=10,I2=1,C1=0,C2=0), parms = NULL, times = seq(0,T,1), func = SI, events = list(func = NULL, root = TRUE, terminalroot = 1), rootfun=rootfun, method = "radau")
    out<- data.frame(out)
    if(max(out$I1)>10 & max(out$C1)<=C1max & max(out$C2)>C2max & max(out$time)<T){
      outcome = "[mit, circ]"
      max_supp_circ2[i] = u2+yinc/2 
    } else if(max(out$I1)>10 & max(out$C1)>C1max & max(out$C2)<=C2max& max(out$time)<T){
      outcome = "[circ, cont]"
    } else if(max(out$I1)>10 & max(out$C1)<=C1max & max(out$C2)<=C2max& max(out$time)<T){
      outcome = "[mit, cont]"
      max_supp2[i] = u2+yinc/2 
    } else if(max(out$time)==T){
      outcome = "no end"
    } else if(max(out$I1)>10 & max(out$C1)>C1max & max(out$C2)>C2max& max(out$time)<T){
      outcome = "[circ, circ]"
      max_circ_mit2[i] = u2+yinc/2 
    }
    else if(max(out$I1)<=10 & max(out$time)<T & max(out$C1)<C1max & max(out$C2)<C2max){
      outcome = "[elim, elim]"
      max_elim2[i] = u2+yinc/2 
    }
    J = cumsum(beta*out$S*(out$I1+c*out$I2)*diff(c(0,out$time)))
    stratMx2 = rbind(stratMx2, data.frame(u1 = u1, u2 = u2, strategy = outcome, J=J, T = max(out$time)))
  }}

elim2 = data.frame(u1vec, max_elim2, max_supp2, max_supp_circ2)

g2=ggplot() +
  geom_raster(data = stratMx2,aes(x=u1, y=u2, fill= strategy),hjust = 1, vjust=1)+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_elim2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp_circ2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_circ_mit2),col="white")+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  ggtitle("High importations")+
  scale_fill_manual(values=cols, breaks = c("[elim, elim]", "[mit, cont]", "[mit, circ]", "[circ, cont]", "[circ, circ]"))+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size

g3=ggplot() +
  geom_raster(data = stratMx,aes(x=u1, y=u2, fill=T),hjust = 1, vjust=1)+
  geom_line(data=elim,aes(x=u1vec+1.5*xinc,y=max_elim-2*yinc),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_supp),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_supp_circ),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_circ_mit),col="white")+
  annotate("text", x = c(1.5, 0.5, .4, 0.05, 0.05), y = c(1.5,1.5,.1,1, 2.7), label = c("[elim, elim]","[circ,circ]", "[circ, cont]", "[mit, cont]", "[mit, circ]"),angle = c(0,0,0,90,90), col="black", size = 5)+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  ggtitle("Low importations")+scale_fill_gradient(low = cols[3], high = "black")+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size


g4=ggplot() +
  geom_raster(data = stratMx2,aes(x=u1, y=u2, fill=T),hjust = 1, vjust=1)+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_elim2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp_circ2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_circ_mit2),col="white")+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  annotate("text", x = c(1.75, 1, 0.05, 0.05), y = c(2,2,0.2, 2), label = c("[elim, elim]","[circ,circ]", "[mit, cont]", "[mit, circ]"),angle = c(0,0,90,90), col=c("black", "white", "black", "black"), size = c(5,5,4,5))+
  ggtitle("High importations")+scale_fill_gradient(low = cols[3], high = "black")+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size

g5=ggplot() +
  geom_raster(data = stratMx,aes(x=u1, y=u2, fill= J),hjust = 1, vjust=1)+
  geom_line(data=elim,aes(x=u1vec+1.5*xinc,y=max_elim-2*yinc),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_supp),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_supp_circ),col="white")+
  geom_line(data=elim,aes(x=u1vec+xinc,y=max_circ_mit),col="white")+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  annotate("text", x = c(1.5, 0.5, .4, 0.05, 0.05), y = c(1.5,1.5,.1,1, 2.7), label = c("[elim, elim]","[circ,circ]", "[circ, cont]", "[mit, cont]", "[mit, circ]"),angle = c(0,0,0,90,90), col=c("black", "white", "white", "white", "white"), size = 5)+
  ggtitle("Low importations")+
  scale_fill_gradient(low = cols[2], high = "black")+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size

g6=ggplot() +
  geom_raster(data = stratMx2,aes(x=u1, y=u2, fill= J),hjust = 1, vjust=1)+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_elim2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_supp_circ2),col="white")+
  geom_line(data=elim2,aes(x=u1vec+xinc,y=max_circ_mit2),col="white")+
  xlab(expression(paste("community isolation daily max,", u[1][max])))+ylab(expression(paste("traveler isolation daily max,", u[2][max])))+
  ggtitle("High importations")+
  annotate("text", x = c(1.75, 1, 0.05, 0.05), y = c(2,2,0.2, 2), label = c("[elim, elim]","[circ,circ]", "[mit, cont]", "[mit, circ]"),angle = c(0,0,90,90), col=c("black", "white", "white", "white"), size = c(5,5,4,5))+
  scale_fill_gradient(low = cols[2], high = "black")+theme_classic()+
  theme(
  plot.title = element_text(size = 30),       # Title size
  axis.title = element_text(size = 30),       # Axis title size
  axis.text = element_text(size = 30),        # Axis text size
  legend.title = element_text(size = 20),     # Legend title size
  legend.text = element_text(size = 20) )      # Legend text size

g8=(g1+g2)+plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 25))


g9 = (g3+g4)+plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 25))

g10 = (g5+g6)+plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 25))


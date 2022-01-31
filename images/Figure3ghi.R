MMcol <- viridis::plasma(n=5, alpha = 0.66)[2]
INHcol <- viridis::plasma(n=5)[4]

mm_f <- function(S, Vmax, Km){
  R = (Vmax*S)/(Km+S)  # Michaelis-Menten function
  return(R)
}

inh_f <- function(S, Vmax, Km, Ki){
  R = Vmax * S / (Km + S *(1+S/Ki))  # 3 Parameter Model
  return(R)
}

Km_i = 1
Ki_i = 1
nd <- data.frame(S = seq(from = 0, to = 10, length.out = 10001)) ### m20 ###

pdf("images/Figure3ghi.pdf", width = 2.5, height = 4.5)
def.p <- function(){
  plot(y = 1, x=1, ylab = "Reaction rate", xlab = "Substrate concentration",
       xlim = c(0, 10), ylim = c(0,2), type = "n")
  lines(x = nd$S, y = mm_f(nd$S, Vmax = 1, Km = Km_i), col = MMcol) # 2P-MM
  lines(x = nd$S, y = inh_f(nd$S, Vmax = 1, Km = Km_i, Ki = Ki_i), col = INHcol) # 3P-MM
}
par(mfrow = c(3, 1), mar = c(2,4,1,2))
def.p();title(main = "Vmax")
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1*2, Km = Km_i, Ki = Ki_i), col = INHcol, lty = 2)
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1/2, Km = Km_i, Ki = Ki_i), col = INHcol, lty = 2)
lines(x = nd$S, y = mm_f(nd$S, Vmax = 1*2, Km = Km_i), col = MMcol, lty = 2) # 2P
lines(x = nd$S, y = mm_f(nd$S, Vmax = 1/2, Km = Km_i), col = MMcol, lty = 2) # 2P
def.p();title(main = "Km")
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1, Km = Km_i*2, Ki = Ki_i), col = INHcol, lty = 2)
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1, Km = Km_i/2, Ki = Ki_i), col = INHcol, lty = 2)
lines(x = nd$S, y = mm_f(nd$S, Vmax = 1, Km = Km_i*2), col = MMcol, lty = 2) # 2P
lines(x = nd$S, y = mm_f(nd$S, Vmax = 1, Km = Km_i/2), col = MMcol, lty = 2) # 2P
def.p();title(main = "Ki")
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1, Km = Km_i, Ki = Ki_i*2), col = INHcol, lty = 2)
lines(x = nd$S, y = inh_f(nd$S, Vmax = 1, Km = Km_i, Ki = Ki_i/2), col = INHcol, lty = 2)
dev.off()

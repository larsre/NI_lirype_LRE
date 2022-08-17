

library(diagram)

Num <- 2

DiffMat <- matrix(data=0, nrow=Num, ncol=Num)
AA <- as.data.frame(DiffMat)

AA[1,c(1,2)] <- "S[t] * F[t]"
AA[2,c(1,2)] <- "S[t]"
rownames(AA) <- colnames(AA) <- c("juvenile", "adult")

plotmat(AA, 2, lwd=3, arr.len=0.5, self.lwd=c(2,2), curve=matrix(c(0.25, 0.25, 0.75, 0.5), ncol=2),  
          my=-0.25, self.shiftx=matrix(c(-0.1, 0.1, 0.5, 0.5), ncol=2), 
          self.shifty = matrix(c(0.15, 0.15, 0, 0.02), ncol=2))

name <- c(expression(A[1]), expression(A[2]))

diagram::plotmat(A=AA, pos=1, curve=0.75, name=name, lwd=2, relsize=0.87, lcol="grey",
        arr.len=0.3, arr.width=0.05, my=-0.25, self.cex=1.5, cex.txt=0.8, self.lwd=c(1,2),
        box.size=0.035, arr.type="simple", dtext=0.8, box.cex=0.8,
        main="", box.lcol="dark red")





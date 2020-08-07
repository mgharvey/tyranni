setwd('~/Documents/research/tyranni/v2/div_analyses/div_through_time')
getwd()

library(ape)
library(phytools)
library(TreePar)

grid <- 2
start <- 2
end <- 40
m <- ((end-start)/grid)+1

# Modified plotting function
bd.shifts.plot.mod <- function (resall, shifts, timemax = 100, ratemin = -1, ratemax = 1, plotturnover = FALSE, linecol = "blue", plotaxes=TRUE) {
    pick <- function(resall, shifts) {
        resall[[shifts + 1]]
    }
    estimatesall <- sapply(resall, pick, shifts = shifts)
    plot(c(0, -timemax), c(ratemin, ratemax), col = "white", 
        xlab = "time before the present", ylab = "diversification rate")
    for (i in 1:length(estimatesall[1, ])) {
        estimates <- estimatesall[, i]
        rates <- length(estimates)/3
        estimates <- estimates[-1]
        if (rates > 1) {
            time <- estimates[(length(estimates) - rates + 2):length(estimates)]
            time <- sort(c(time, time, 0, timemax))
            turnover <- estimates[1]
            div <- estimates[rates + 1]
            for (j in 1:(rates - 1)) {
                turnover <- c(turnover, estimates[j:(j + 1)])
                div <- c(div, estimates[(rates + j):(rates + 
                  j + 1)])
            }
            turnover <- c(turnover, estimates[rates])
            div <- c(div, estimates[2 * rates])
        }
        else {
            time <- c(0, timemax)
            turnover <- c(estimates[1], estimates[1])
            div <- c(estimates[2], estimates[2])
        }
        if(plotaxes==TRUE) {
        		lines(-time, div, col=linecol, ylab="Time before the present (My)", xlab="Diversification rate")
        } else {
       	 	lines(-time, div, col=linecol, xaxt="n", yaxt="n", ylab="", xlab="")
        }
    }
}

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
x <- sort(getx(tree), decreasing=TRUE)
res1 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res1

# Likelihood ratio test
for(i in 1:m) {
	print(i)
	test <- pchisq(2*(res1[[2]][[i]][1]-res1[[2]][[i+1]][1]),3)
	print(test) # check if >0.95
	if(test > 0.95) {
		best <- i
	}
}
print(best) # best i (= # shifts-1)
res[[2]][[best]] # parameters for that model

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
x <- sort(getx(tree), decreasing=TRUE)
res2 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res2

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
x <- sort(getx(tree), decreasing=TRUE)
res3 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res3

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
x <- sort(getx(tree), decreasing=TRUE)
res4 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res4

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
x <- sort(getx(tree), decreasing=TRUE)
res5 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res5

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
x <- sort(getx(tree), decreasing=TRUE)
res6 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res6

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
x <- sort(getx(tree), decreasing=TRUE)
res7 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res7

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
x <- sort(getx(tree), decreasing=TRUE)
res8 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res8

tree <- read.tree('~/Documents/research/tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
x <- sort(getx(tree), decreasing=TRUE)
res9 <- bd.shifts.optim(x, sampling=c(rep(1, m+1)), grid, start, end, yule=FALSE)
res9

# Plot

pdf("divTTplots.pdf", width=6, height=6)

bd.shifts.plot.mod(list(res1[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="black", plotaxes=TRUE)
par(new=TRUE)
bd.shifts.plot.mod(list(res2[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="purple", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res3[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="blue", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res4[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="green", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res5[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="red", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res6[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="orange", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res7[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="yellow", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res8[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="gray", plotaxes=FALSE)
par(new=TRUE)
bd.shifts.plot.mod(list(res9[[2]]), m-1, timemax=max(nodeHeights(tree)), ratemin=-.4, ratemax=1.1, linecol="lightblue", plotaxes=FALSE)


legend("topleft", c("T400F Howard & Moore (H&M)", "T400F Clements", "T400F IUCN", "T400F 1 My cutoff", "T400F 2 My cutoff", "HGAPF H&M", "T400F Astral H&M", "T400F Astral H&M poor samples dropped"), lty=1, col=c("black", 
"purple", "blue", "green", "red", "orange", "yellow", "gray", "lightblue"), bty="n")

dev.off()


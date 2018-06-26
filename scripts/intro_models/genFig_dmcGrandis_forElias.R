cl_mig_stagSweeps = readRDS("~/analysis/data/intro_models/compLikelihood_mig_stagSweeps.RDS")
cl_neutral = readRDS("~/analysis/data/intro_models/compLikelihood_neutral.RDS")
cl_sv = readRDS("~/analysis/data/intro_models/compLikelihood_stdVar.RDS")

sels = readRDS("~/analysis/data/intro_models/sels.RDS")
times = readRDS("~/analysis/data/intro_models/times.RDS")
gs = readRDS("~/analysis/data/intro_models/gs.RDS")


maxCL_mig_stagSweeps.sels = sapply(1 : length(sels), function(i) max(unlist(cl_mig_stagSweeps[[i]])))
maxCL_sv.sels = sapply(1 : length(sels), function(i) max(unlist(cl_sv[[i]])))

cl_mig_stagSweeps_byTime = lapply(1 : length(times), function(time) lapply(1 : length(sels), function(sel) lapply(1 : length(gs), function(g) cl_mig_stagSweeps[[sel]][[g]][[time]])))
maxCL_mig_stagSweeps.times = sapply(1 : length(times), function(i) max(unlist(cl_mig_stagSweeps_byTime[[i]])))

sels18_scaled = (sels[18:31] - sels[18]) / (sels[31] - sels[18]) #get last sels and put on (0,1) scale
times18_scaled = (times[1:18] - times[1]) / (times[18] - times[1]) #get last times and put on (0,1) scale

yMin = min(c(maxCL_mig_stagSweeps.sels, maxCL_sv.sels) - cl_neutral)
yMax = max(c(maxCL_mig_stagSweeps.sels, maxCL_sv.sels) - cl_neutral)

pdf("compLike_grandis_sels.pdf", width = 5, height = 5)
par(oma = c(1,1.5,2,1.5))
plot(sels18_scaled, maxCL_mig_stagSweeps.sels[18:31] - cl_neutral, cex.axis = 1.3, ylab = "Composite log-likelihood (model - neutral)", ylim = c(4.8e5, yMax), type = "l", col = "red", xlim = c(0, 1), xaxt = "none", xlab = "", lwd = 3, cex.lab = 1.1)
lines(sels18_scaled, maxCL_sv.sels[18:31] - cl_neutral, col = "red", lwd = 3, lty = 2)
axis(1, at = seq(0, 1, length.out = 6), labels = seq(0.5, 1, length.out = 6), col = "red", lwd = 2, cex.axis = 1.3)
lines(times18_scaled, maxCL_mig_stagSweeps.times[1:18] - cl_neutral, col = "plum2", lwd = 3, lty = 1)
axis(3, lwd = 2, at = seq(0, 1, length.out = 6), labels = seq(times[1], times[18], length.out = 6), col = "plum2", cex.axis = 1.3)
legend("bottomright", lty = c(1,2), legend = c("Introgression", "ILS"), lwd = 3)
mtext(line = 3, side = 1, text = "selection coefficient", cex = 1.3)
mtext(line = 3, side = 3, text = "time between introgression and selection (generations)", cex = 1.1)
dev.off() 
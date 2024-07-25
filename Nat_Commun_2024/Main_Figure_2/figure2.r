library(corrplot)

####################
corrdf <- data.frame(
	'MCF7'=c(0, 0.65, 0.67, 0.43, 0.20, 0.39),
	'MCF7M1'=c(0.65, 0, 0.86, 0.25, 0.61, 0.34),
	'MCF7TR'=c(0.67, 0.86, 0, 0.13, 0.41, 0.58),
	'scMCF7'=c(0.43, 0.25, 0.13, 0, 0.05, 0.07),
	'scMCF7M1'=c(0.20, 0.61, 0.41, 0.05, 0, 0.28),
	'scMCF7TR'=c(0.39, 0.34, 0.58, 0.07, 0.28, 0))
rownames(corrdf) <- colnames(corrdf)
abc <- as.matrix(corrdf)
#corrplot(abc, order = "AOE", type = "upper", tl.pos = "d", method = "pie")
#corrplot(abc, add = TRUE, type = "lower", method = "number", order = "AOE",
#         diag = FALSE, tl.pos = "n", cl.pos = "n")
corrplot(abc, type = "upper", tl.pos = "d", method = "pie", cl.lim=c(0,1))
corrplot(abc, add = TRUE, type = "lower", method = "number", 
         diag = FALSE, tl.pos = "n", cl.pos = "n")

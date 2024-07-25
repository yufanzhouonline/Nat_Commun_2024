### inhibitor figures

########################################################################
### inhibitor figures with error bars
### Curve of Real time live cell counts

### LDN-192960
ldnmcf7 <- c(100, 66.44910583, 60.30946974, 46.53029271, 43.26082054, 50.88551728, 62.6259015)
ldnmcf7m1 <- c(100, 69.65389043, 64.40470993, 46.77555513, 37.41148539, 34.79123934, 39.60995317)
ldnmcf7tr <- c(100, 60.32417118, 38.79595905, 22.42126141, 13.13265712, 4.844496187, 3.714643094)
ebarmcf7 <- c(0, 1.715576517, 1.20081294, 1.580952226, 2.607470565, 3.091631769, 4.254833465)
ebarmcf7m1 <- c(0, 0.945199681, 1.29347574, 1.937080228, 1.688257085, 2.25025786, 3.305552312)
ebarmcf7tr <- c(0, 1.663625374, 3.159788585, 1.849223023, 1.305803321, 0.34313094, 0.322233416)
par(mar = c(5, 10, 4, 2) + 0.1)
plot(0:6, ldnmcf7, type="o", col="green", xlim=c(0, 6.5), ylim=c(0,110), lwd=5, cex.lab=2, cex.axis=2, xlab="Day", ylab="Real-time live cell counts (%)")
lines(0:6, ldnmcf7m1, type="o", col="blue", lwd=5)
lines(0:6, ldnmcf7tr, type="o", col="red", lwd=5)
for (i in 1:6)
{
    arrows(i, ldnmcf7[i+1], i, ldnmcf7[i+1] - ebarmcf7[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7[i+1], i, ldnmcf7[i+1] + ebarmcf7[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7m1[i+1], i, ldnmcf7m1[i+1] - ebarmcf7m1[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7m1[i+1], i, ldnmcf7m1[i+1] + ebarmcf7m1[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7tr[i+1], i, ldnmcf7tr[i+1] - ebarmcf7tr[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7tr[i+1], i, ldnmcf7tr[i+1] + ebarmcf7tr[i+1], angle=90, length=0.1)
}
lines(c(6.1, 6.3), c(ldnmcf7[7], ldnmcf7[7]), lwd=5)
lines(c(6.1, 6.3), c(ldnmcf7tr[7], ldnmcf7tr[7]), lwd=5)
lines(c(6.3, 6.3), c(ldnmcf7[7], ldnmcf7tr[7]), lwd=5)
text(6.4, (ldnmcf7[7] + ldnmcf7tr[7]) / 2, "*", cex=2)
text(6.5, (ldnmcf7[7] + ldnmcf7tr[7]) / 2, "*", cex=2)
legend(4.3, 110, legend=c("MCF7", "MCF7M1", "MCF7TR"), col=c("green", "blue", "red"), lty=1, cex=2, lwd=5)

t.test(ldnmcf7[2:7], ldnmcf7m1[2:7], paired = TRUE, alternative = "two.sided")
### p-value = 0.2275
t.test(ldnmcf7[2:7], ldnmcf7tr[2:7], paired = TRUE, alternative = "two.sided")
### p-value = 0.009702

### MS023
ldnmcf7 <- c(100, 93.27904128, 87.47567732, 78.6703252, 67.9733854, 62.47543204, 62.41508412)
ldnmcf7m1 <- c(100, 91.2585448, 83.45378705, 73.19969463, 68.80942386, 67.56309873, 65.71305791)
ldnmcf7tr <- c(100, 89.98438495, 82.44646953, 69.00512459, 61.95679312, 57.98897828, 41.1966133)
ebarmcf7 <- c(0, 1.9408374,2.455932649, 3.28387141, 4.997263747, 5.547125727, 6.20506964)
ebarmcf7m1 <- c(0, 2.069573762, 2.559664811, 4.3352897, 3.055534838, 2.87540977, 3.616378266)
ebarmcf7tr <- c(0, 1.936980535, 1.937962772, 3.741928302, 4.109714392, 5.219456446, 4.044205307)
par(mar = c(5, 10, 4, 2) + 0.1)
plot(0:6, ldnmcf7, type="o", col="green", xlim=c(0, 6.5), ylim=c(0,119), lwd=5, cex.lab=2, cex.axis=2, xlab="Day", ylab="Real-time live cell counts (%)")
lines(0:6, ldnmcf7m1, type="o", col="blue", lwd=5)
lines(0:6, ldnmcf7tr, type="o", col="red", lwd=5)
for (i in 1:6)
{
    arrows(i, ldnmcf7[i+1], i, ldnmcf7[i+1] - ebarmcf7[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7[i+1], i, ldnmcf7[i+1] + ebarmcf7[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7m1[i+1], i, ldnmcf7m1[i+1] - ebarmcf7m1[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7m1[i+1], i, ldnmcf7m1[i+1] + ebarmcf7m1[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7tr[i+1], i, ldnmcf7tr[i+1] - ebarmcf7tr[i+1], angle=90, length=0.1)
    arrows(i, ldnmcf7tr[i+1], i, ldnmcf7tr[i+1] + ebarmcf7tr[i+1], angle=90, length=0.1)
}
lines(c(6.1, 6.3), c(ldnmcf7[7], ldnmcf7[7]), lwd=5)
lines(c(6.1, 6.3), c(ldnmcf7tr[7], ldnmcf7tr[7]), lwd=5)
lines(c(6.3, 6.3), c(ldnmcf7[7], ldnmcf7tr[7]), lwd=5)
text(6.4, (ldnmcf7[7] + ldnmcf7tr[7]) / 2, "*", cex=2)
legend(4, 119, legend=c("MCF7", "MCF7M1", "MCF7TR"), col=c("green", "blue", "red"), lty=1, cex=2, lwd=5)

t.test(ldnmcf7[2:7], ldnmcf7m1[2:7], paired = TRUE, alternative = "two.sided")
### p-value = 0.8312
t.test(ldnmcf7[2:7], ldnmcf7tr[2:7], paired = TRUE, alternative = "two.sided")
### p-value = 0.02907





########################################################################
### cell proliferation bars


### MCF7
dmsomcf7 <- c(100, 135.9922179, 271.4007782, 253.1128405)
ms023mcf7 <- c(100, 124.8262548, 249.1891892, 232.0849421)
ldnmcf7 <- c(100, 115.7692308, 221.6923077, 215.6538462)
error1 <- c(6.137384564, 16.72441257, 14.7594991, 44.91454032)
error2 <- c(8.240507535, 5.057325557, 23.39147056, 35.323881)
error3 <- c(8.760370017, 5.636551256, 23.08451798, 27.52908629)

mcf7 <- t(data.frame(dmsomcf7, ms023mcf7, ldnmcf7))
colnames(mcf7) <- c("0", "2", "4", "6")
cells <- as.matrix(mcf7)
par(mar = c(5, 10, 4, 2) + 0.1)
barplot(cells, beside=TRUE, col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), ylim=c(0, 300), xlab="Day", ylab="MCF7 cell proliferation (%)", cex.names=2, cex.axis=2, cex.lab=2)
legend(1, 300, legend=c("DMSO", "MS023", "LDN-192960"), col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), pch=15, cex=2)


dmsomcf7d0 <- c(92.91828794, 103.7743191, 103.307393)
ms023mcf7d0 <- c(90.57915058, 105.8687259, 103.5521236)
ldnmcf7d0 <- c(89.88461538, 105, 105.1153846)

dmsomcf7d2 <- c(123.618677, 129.3385214, 155.0194553)
ms023mcf7d2 <- c(122.2007722, 121.6216216, 130.6563707)
ldnmcf7d2 <- c(119.3076923, 118.7307692, 109.2692308)

dmsomcf7d4 <- c(254.7081712, 282.7237354, 276.770428)
ms023mcf7d4 <- c(222.5096525, 258.8803089, 266.1776062)
ldnmcf7d4 <- c(209.6538462, 207.1153846, 248.3076923)

dmsomcf7d6 <- c(202.1789883, 270.1167315, 287.0428016)
ms023mcf7d6 <- c(192.2779923, 244.2857143, 259.6911197)
ldnmcf7d6 <- c(184.6153846, 225.2307692, 237.1153846)

points(c(1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5), dmsomcf7d0, col = "black", pch = 20)
points(c(2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5), ms023mcf7d0, col = "black", pch = 20)
points(c(3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5), ldnmcf7d0, col = "black", pch = 20)

points(c(5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5), dmsomcf7d2, col = "black", pch = 20)
points(c(6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5), ms023mcf7d2, col = "black", pch = 20)
points(c(7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5), ldnmcf7d2, col = "black", pch = 20)

points(c(9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5), dmsomcf7d4, col = "black", pch = 20)
points(c(10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5), ms023mcf7d4, col = "black", pch = 20)
points(c(11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5), ldnmcf7d4, col = "black", pch = 20)

points(c(13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5), dmsomcf7d6, col = "black", pch = 20)
points(c(14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5), ms023mcf7d6, col = "black", pch = 20)
points(c(15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5), ldnmcf7d6, col = "black", pch = 20)

for (i in 1:4)
{
    arrows(i * 4 - 2.5, dmsomcf7[i], i * 4 - 2.5, dmsomcf7[i] + error1[i], angle=90, length=0.1)
    arrows(i * 4 - 1.5, ms023mcf7[i], i * 4 - 1.5, ms023mcf7[i] + error2[i], angle=90, length=0.1)
    arrows(i * 4 - 0.5, ldnmcf7[i], i * 4 - 0.5, ldnmcf7[i] + error3[i], angle=90, length=0.1)
}

### MCF7M1
dmsomcf7 <- c(100, 127.4611399, 241.7696293, 236.6281387)
ms023mcf7 <- c(100, 121.2688085, 208.6620577, 223.7088247)
ldnmcf7 <- c(100, 130.2364865, 228.8006757, 238.597973)
error1 <- c(7.509437348, 10.97697959, 12.7670943, 22.59231756)
error2 <- c(6.466426746, 2.547451771, 13.96061087, 10.56158817)
error3 <- c(5.459897419, 13.05636435, 20.84297071, 28.11082474)

mcf7 <- t(data.frame(dmsomcf7, ms023mcf7, ldnmcf7))
colnames(mcf7) <- c("0", "2", "4", "6")
cells <- as.matrix(mcf7)
par(mar = c(5, 10, 4, 2) + 0.1)
barplot(cells, beside=TRUE, col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), ylim=c(0, 300), xlab="Day", ylab="MCF7M1 cell proliferation (%)", cex.names=2, cex.axis=2, cex.lab=2)
legend(1, 300, legend=c("DMSO", "MS023", "LDN-192960"), col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), pch=15, cex=2)


dmsomcf7d0 <- c(91.35113591, 104.862495, 103.7863691)
ms023mcf7d0 <- c(92.59861732, 104.554697, 102.8466856)
ldnmcf7d0 <- c(95.77702703, 106.1655405, 98.05743243)

dmsomcf7d2 <- c(139.657234, 118.3738541, 124.3523316)
ms023mcf7d2 <- c(120.0488003, 119.5607971, 124.196828)
ldnmcf7d2 <- c(145.3125, 122.7618243, 122.6351351)

dmsomcf7d4 <- c(253.2483061, 244.0414508, 228.0191311)
ms023mcf7d4 <- c(224.6034974, 198.6173241, 202.7653518)
ldnmcf7d4 <- c(249.8310811, 208.1503378, 228.4206081)

dmsomcf7d6 <- c(214.9860502, 260.0637704, 234.8345955)
ms023mcf7d6 <- c(214.2334282, 221.7974786, 235.0955673)
ldnmcf7d6 <- c(270.9881757, 220.5658784, 224.2398649)

points(c(1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5), dmsomcf7d0, col = "black", pch = 20)
points(c(2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5), ms023mcf7d0, col = "black", pch = 20)
points(c(3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5), ldnmcf7d0, col = "black", pch = 20)

points(c(5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5), dmsomcf7d2, col = "black", pch = 20)
points(c(6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5), ms023mcf7d2, col = "black", pch = 20)
points(c(7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5), ldnmcf7d2, col = "black", pch = 20)

points(c(9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5), dmsomcf7d4, col = "black", pch = 20)
points(c(10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5), ms023mcf7d4, col = "black", pch = 20)
points(c(11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5), ldnmcf7d4, col = "black", pch = 20)

points(c(13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5), dmsomcf7d6, col = "black", pch = 20)
points(c(14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5), ms023mcf7d6, col = "black", pch = 20)
points(c(15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5), ldnmcf7d6, col = "black", pch = 20)

for (i in 1:4)
{
    arrows(i * 4 - 2.5, dmsomcf7[i], i * 4 - 2.5, dmsomcf7[i] + error1[i], angle=90, length=0.1)
    arrows(i * 4 - 1.5, ms023mcf7[i], i * 4 - 1.5, ms023mcf7[i] + error2[i], angle=90, length=0.1)
    arrows(i * 4 - 0.5, ldnmcf7[i], i * 4 - 0.5, ldnmcf7[i] + error3[i], angle=90, length=0.1)
}

for (i in 1:4)
{
    arrows(i * 4 - 2.5, dmsomcf7[i], i * 4 - 2.5, dmsomcf7[i] + error1[i], angle=90, length=0.1)
    arrows(i * 4 - 1.5, ms023mcf7[i], i * 4 - 1.5, ms023mcf7[i] + error2[i], angle=90, length=0.1)
    arrows(i * 4 - 0.5, ldnmcf7[i], i * 4 - 0.5, ldnmcf7[i] + error3[i], angle=90, length=0.1)
}

### MCF7TR
dmsomcf7 <- c(100, 101.7746914, 128.0092593, 166.0493827)
ms023mcf7 <- c(100, 94.54545455, 103.3333333, 93.03030303)
ldnmcf7 <- c(100, 92.55236618, 87.58727696, 86.42358417)
error1 <- c(2.891202777, 4.677606348, 13.57476566, 25.52741562)
error2 <- c(6.661068586, 3.363329224, 12.6383528, 19.63816981)
error3 <- c(5.50923841, 1.807785911, 4.905946239, 1.148072039)

mcf7 <- t(data.frame(dmsomcf7, ms023mcf7, ldnmcf7))
colnames(mcf7) <- c("0", "2", "4", "6")
cells <- as.matrix(mcf7)
par(mar = c(5, 10, 4, 2) + 0.1)
barplot(cells, beside=TRUE, col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), ylim=c(0, 200), xlab="Day", ylab="MCF7TR cell proliferation (%)", cex.names=2, cex.axis=2, cex.lab=2)
legend(1, 200, legend=c("DMSO", "MS023", "LDN-192960"), col=c("#88FF88AA", "#8888FFAA", "#FF8888AA"), pch=15, cex=2)


dmsomcf7d0 <- c(103.2407407, 99.07407407, 97.68518519)
ms023mcf7d0 <- c(107.5, 94.77272727, 97.72727273)
ldnmcf7d0 <- c(106.3615206, 96.81923972, 96.81923972)

dmsomcf7d2 <- c(99.07407407, 99.07407407, 107.1759259)
ms023mcf7d2 <- c(92.27272727, 92.95454545, 98.40909091)
ldnmcf7d2 <- c(94.02637704, 93.09542281, 90.53529868)

dmsomcf7d4 <- c(118.287037, 122.2222222, 143.5185185)
ms023mcf7d4 <- c(89.31818182, 106.8181818, 113.8636364)
ldnmcf7d4 <- c(82.15671063, 91.69899147, 88.90612878)

dmsomcf7d6 <- c(136.5740741, 180.5555556, 181.0185185)
ms023mcf7d6 <- c(70.90909091, 99.77272727, 108.4090909)
ldnmcf7d6 <- c(85.64778898, 87.742436, 85.88052754)

points(c(1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5, 1.25 + runif(1) * 0.5), dmsomcf7d0, col = "black", pch = 20)
points(c(2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5, 2.25 + runif(1) * 0.5), ms023mcf7d0, col = "black", pch = 20)
points(c(3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5, 3.25 + runif(1) * 0.5), ldnmcf7d0, col = "black", pch = 20)

points(c(5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5, 5.25 + runif(1) * 0.5), dmsomcf7d2, col = "black", pch = 20)
points(c(6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5, 6.25 + runif(1) * 0.5), ms023mcf7d2, col = "black", pch = 20)
points(c(7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5, 7.25 + runif(1) * 0.5), ldnmcf7d2, col = "black", pch = 20)

points(c(9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5, 9.25 + runif(1) * 0.5), dmsomcf7d4, col = "black", pch = 20)
points(c(10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5, 10.25 + runif(1) * 0.5), ms023mcf7d4, col = "black", pch = 20)
points(c(11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5, 11.25 + runif(1) * 0.5), ldnmcf7d4, col = "black", pch = 20)

points(c(13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5, 13.25 + runif(1) * 0.5), dmsomcf7d6, col = "black", pch = 20)
points(c(14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5, 14.25 + runif(1) * 0.5), ms023mcf7d6, col = "black", pch = 20)
points(c(15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5, 15.25 + runif(1) * 0.5), ldnmcf7d6, col = "black", pch = 20)


for (i in 1:4)
{
    arrows(i * 4 - 2.5, dmsomcf7[i], i * 4 - 2.5, dmsomcf7[i] + error1[i], angle=90, length=0.1)
    arrows(i * 4 - 1.5, ms023mcf7[i], i * 4 - 1.5, ms023mcf7[i] + error2[i], angle=90, length=0.1)
    arrows(i * 4 - 0.5, ldnmcf7[i], i * 4 - 0.5, ldnmcf7[i] + error3[i], angle=90, length=0.1)
}




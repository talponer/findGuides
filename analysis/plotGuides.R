######################################################################
### Script to plot the guides location in the CDS
######################################################################

## AARS
AARSguides <- read.table("../AARS_putativeGuides.log", header=T, as.is=T, sep="\t")
AARS.l <- 2907
AARSguides.good <- AARSguides[which(AARSguides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- AARSguides.good$To[which(AARSguides.good$Strand == "+")[1]]
for (I in which(AARSguides.good$Strand == "+")[-1]){
    startPos <- AARSguides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- AARSguides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- AARSguides.good$From[which(AARSguides.good$Strand == "-")[1]]
for (I in which(AARSguides.good$Strand == "-")[-1]){
    startPos <- AARSguides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- AARSguides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/AARSguidesDist.pdf", width=20)
plot(c(1,AARS.l), c(1,1), type="l", xlim=c(1,AARS.l), ylim=c(0,2), main="AARS", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(AARSguides.good$From[which(AARSguides.good$Strand == "+")], yposPlus, AARSguides.good$To[which(AARSguides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(AARSguides.good$From[which(AARSguides.good$Strand == "-")], yposMinus, AARSguides.good$To[which(AARSguides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()

## DAD1
DAD1guides <- read.table("../DAD1_putativeGuides.log", header=T, as.is=T, sep="\t")
DAD1.l <- 342
DAD1guides.good <- DAD1guides[which(DAD1guides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- DAD1guides.good$To[which(DAD1guides.good$Strand == "+")[1]]
for (I in which(DAD1guides.good$Strand == "+")[-1]){
    startPos <- DAD1guides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- DAD1guides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- DAD1guides.good$From[which(DAD1guides.good$Strand == "-")[1]]
for (I in which(DAD1guides.good$Strand == "-")[-1]){
    startPos <- DAD1guides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- DAD1guides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/DAD1guidesDist.pdf", width=20)
plot(c(1,DAD1.l), c(1,1), type="l", xlim=c(1,DAD1.l), ylim=c(0,2), main="DAD1", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(DAD1guides.good$From[which(DAD1guides.good$Strand == "+")], yposPlus, DAD1guides.good$To[which(DAD1guides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(DAD1guides.good$From[which(DAD1guides.good$Strand == "-")], yposMinus, DAD1guides.good$To[which(DAD1guides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()


## DDX3X
DDX3Xguides <- read.table("../DDX3X_putativeGuides.log", header=T, as.is=T, sep="\t")
DDX3X.l <- 1989
DDX3Xguides.good <- DDX3Xguides[which(DDX3Xguides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- DDX3Xguides.good$To[which(DDX3Xguides.good$Strand == "+")[1]]
for (I in which(DDX3Xguides.good$Strand == "+")[-1]){
    startPos <- DDX3Xguides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- DDX3Xguides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- DDX3Xguides.good$From[which(DDX3Xguides.good$Strand == "-")[1]]
for (I in which(DDX3Xguides.good$Strand == "-")[-1]){
    startPos <- DDX3Xguides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- DDX3Xguides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/DDX3XguidesDist.pdf", width=20)
plot(c(1,DDX3X.l), c(1,1), type="l", xlim=c(1,DDX3X.l), ylim=c(0,2), main="DDX3X", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(DDX3Xguides.good$From[which(DDX3Xguides.good$Strand == "+")], yposPlus, DDX3Xguides.good$To[which(DDX3Xguides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(DDX3Xguides.good$From[which(DDX3Xguides.good$Strand == "-")], yposMinus, DDX3Xguides.good$To[which(DDX3Xguides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()


## HARS
HARSguides <- read.table("../HARS_putativeGuides.log", header=T, as.is=T, sep="\t")
HARS.l <- 1530
HARSguides.good <- HARSguides[which(HARSguides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- HARSguides.good$To[which(HARSguides.good$Strand == "+")[1]]
for (I in which(HARSguides.good$Strand == "+")[-1]){
    startPos <- HARSguides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- HARSguides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- HARSguides.good$From[which(HARSguides.good$Strand == "-")[1]]
for (I in which(HARSguides.good$Strand == "-")[-1]){
    startPos <- HARSguides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- HARSguides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/HARSguidesDist.pdf", width=20)
plot(c(1,HARS.l), c(1,1), type="l", xlim=c(1,HARS.l), ylim=c(0,2), main="HARS", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(HARSguides.good$From[which(HARSguides.good$Strand == "+")], yposPlus, HARSguides.good$To[which(HARSguides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(HARSguides.good$From[which(HARSguides.good$Strand == "-")], yposMinus, HARSguides.good$To[which(HARSguides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()


## HCFC1
HCFC1guides <- read.table("../HCFC1_putativeGuides.log", header=T, as.is=T, sep="\t")
HCFC1.l <- 6108
HCFC1guides.good <- HCFC1guides[which(HCFC1guides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- HCFC1guides.good$To[which(HCFC1guides.good$Strand == "+")[1]]
for (I in which(HCFC1guides.good$Strand == "+")[-1]){
    startPos <- HCFC1guides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- HCFC1guides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- HCFC1guides.good$From[which(HCFC1guides.good$Strand == "-")[1]]
for (I in which(HCFC1guides.good$Strand == "-")[-1]){
    startPos <- HCFC1guides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- HCFC1guides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/HCFC1guidesDist.pdf", width=20)
plot(c(1,HCFC1.l), c(1,1), type="l", xlim=c(1,HCFC1.l), ylim=c(0,2), main="HCFC1", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(HCFC1guides.good$From[which(HCFC1guides.good$Strand == "+")], yposPlus, HCFC1guides.good$To[which(HCFC1guides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(HCFC1guides.good$From[which(HCFC1guides.good$Strand == "-")], yposMinus, HCFC1guides.good$To[which(HCFC1guides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()



## KARS
KARSguides <- read.table("../KARS_putativeGuides.log", header=T, as.is=T, sep="\t")
KARS.l <- 1794
KARSguides.good <- KARSguides[which(KARSguides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- KARSguides.good$To[which(KARSguides.good$Strand == "+")[1]]
for (I in which(KARSguides.good$Strand == "+")[-1]){
    startPos <- KARSguides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- KARSguides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- KARSguides.good$From[which(KARSguides.good$Strand == "-")[1]]
for (I in which(KARSguides.good$Strand == "-")[-1]){
    startPos <- KARSguides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- KARSguides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/KARSguidesDist.pdf", width=20)
plot(c(1,KARS.l), c(1,1), type="l", xlim=c(1,KARS.l), ylim=c(0,2), main="KARS", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(KARSguides.good$From[which(KARSguides.good$Strand == "+")], yposPlus, KARSguides.good$To[which(KARSguides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(KARSguides.good$From[which(KARSguides.good$Strand == "-")], yposMinus, KARSguides.good$To[which(KARSguides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()


## RCC1
RCC1guides <- read.table("../RCC1_putativeGuides.log", header=T, as.is=T, sep="\t")
RCC1.l <- 1266
RCC1guides.good <- RCC1guides[which(RCC1guides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- RCC1guides.good$To[which(RCC1guides.good$Strand == "+")[1]]
for (I in which(RCC1guides.good$Strand == "+")[-1]){
    startPos <- RCC1guides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- RCC1guides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- RCC1guides.good$From[which(RCC1guides.good$Strand == "-")[1]]
for (I in which(RCC1guides.good$Strand == "-")[-1]){
    startPos <- RCC1guides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- RCC1guides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/RCC1guidesDist.pdf", width=20)
plot(c(1,RCC1.l), c(1,1), type="l", xlim=c(1,RCC1.l), ylim=c(0,2), main="RCC1", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(RCC1guides.good$From[which(RCC1guides.good$Strand == "+")], yposPlus, RCC1guides.good$To[which(RCC1guides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(RCC1guides.good$From[which(RCC1guides.good$Strand == "-")], yposMinus, RCC1guides.good$To[which(RCC1guides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()



## RPS4X
RPS4Xguides <- read.table("../RPS4X_putativeGuides.log", header=T, as.is=T, sep="\t")
RPS4X.l <- 792
RPS4Xguides.good <- RPS4Xguides[which(RPS4Xguides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- RPS4Xguides.good$To[which(RPS4Xguides.good$Strand == "+")[1]]
for (I in which(RPS4Xguides.good$Strand == "+")[-1]){
    startPos <- RPS4Xguides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- RPS4Xguides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- RPS4Xguides.good$From[which(RPS4Xguides.good$Strand == "-")[1]]
for (I in which(RPS4Xguides.good$Strand == "-")[-1]){
    startPos <- RPS4Xguides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- RPS4Xguides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/RPS4XguidesDist.pdf", width=20)
plot(c(1,RPS4X.l), c(1,1), type="l", xlim=c(1,RPS4X.l), ylim=c(0,2), main="RPS4X", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(RPS4Xguides.good$From[which(RPS4Xguides.good$Strand == "+")], yposPlus, RPS4Xguides.good$To[which(RPS4Xguides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(RPS4Xguides.good$From[which(RPS4Xguides.good$Strand == "-")], yposMinus, RPS4Xguides.good$To[which(RPS4Xguides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()


## TAF1
TAF1guides <- read.table("../TAF1_putativeGuides.log", header=T, as.is=T, sep="\t")
TAF1.l <- 5619
TAF1guides.good <- TAF1guides[which(TAF1guides$Print == 1),]
yposPlus <- 1.01
newY <- yposPlus
endPos <- TAF1guides.good$To[which(TAF1guides.good$Strand == "+")[1]]
for (I in which(TAF1guides.good$Strand == "+")[-1]){
    startPos <- TAF1guides.good$From[I]
    if(startPos < endPos+5){
        newY <- newY + 0.05
#        yposPlus <- cbind(yposPlus, newY)
    }else{
        newY <- 1.01
        endPos <- TAF1guides.good$To[I]
    }
    yposPlus <- cbind(yposPlus, newY)
}
###
yposMinus <- 0.99
newY <- yposMinus
endPos <- TAF1guides.good$From[which(TAF1guides.good$Strand == "-")[1]]
for (I in which(TAF1guides.good$Strand == "-")[-1]){
    startPos <- TAF1guides.good$To[I]
    if(startPos < endPos+5){
        newY <- newY - 0.05
#        yposMinus <- cbind(yposMinus, newY)
    }else{
        newY <- 0.99
        endPos <- TAF1guides.good$From[I]
    }
    yposMinus <- cbind(yposMinus, newY)
}
pdf("figures/2019-02-12/TAF1guidesDist.pdf", width=20)
plot(c(1,TAF1.l), c(1,1), type="l", xlim=c(1,TAF1.l), ylim=c(0,2), main="TAF1", frame.plot="F", yaxt="n", ylab="", xlab="")
rect(TAF1guides.good$From[which(TAF1guides.good$Strand == "+")], yposPlus, TAF1guides.good$To[which(TAF1guides.good$Strand == "+")], yposPlus+0.04, col="red")
rect(TAF1guides.good$From[which(TAF1guides.good$Strand == "-")], yposMinus, TAF1guides.good$To[which(TAF1guides.good$Strand == "-")], yposMinus-0.04, col="blue")
dev.off()

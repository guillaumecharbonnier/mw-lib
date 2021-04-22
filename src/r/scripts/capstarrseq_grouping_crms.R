dat <- read.table(
    snakemake@input[["fc"]],
    stringsAsFactors = FALSE
)

print("CALCUL DU SEUIL EN FONCTION DU FDR")
for (
    fdr in c(
        as.numeric(snakemake@params[["th_negative_1"]]),
        as.numeric(snakemake@params[["th_negative_2"]])
    )
) {
    print(paste("Compute threshold for fdr=",fdr))
    idx <- which(dat[,4] == "Negative")
    print("Skipping FDR computation if no Negative set is provided.")
    if (length(idx)!=0) {
        P <- ecdf(dat[idx,5])
        th_FoldChange <- quantile(
            P,
            probs = 1-fdr
        )

        print("Identify active and inactive regions")
        idx <- which(dat[,5] >= th_FoldChange)
        groups <- rep(
            "Inactive",
            nrow(dat)
        )
        groups[idx] <- "Active"
        groups[which(is.na(dat[,5]))] <- "NA"
        
        pdf(
            paste0(
                snakemake@params[["outdir"]],
                "Groups.FDR=",
                fdr,
                ".pdf"
            )
        )
        par(mfrow=c(2,2))
        
        print("Produce boxplot for categories")
        categories <- unique(dat[,4])
        fc <- list() ; lfc <- list()
        for (category in categories) {
            idx <- which(dat[,4] == category)
            fc[[category]] <- dat[idx,5]
            lfc[[category]] <- log2(dat[idx,5])
        }
        colors <- rep(
            "lightgrey",
            length(categories)
        )
        if ("Random" %in% categories) {
            colors[which(categories == "Random")] <- "darkblue"
        }
        if ("PosEpromoter" %in% categories) {
            colors[which(categories == "PosEpromoter")] <- "darkgreen"
        }
        if ("Positive" %in% categories) {
            colors[which(categories == "Positive")] <- "darkgreen"
        }
        if ("Negative" %in% categories) {
            colors[which(categories == "Negative")] <- "darkred"
        }
        boxplot(
            fc,
            pch = 20,
            col=colors,
            main = "Fold Change of categories",
            ylab = "Fold Change",
            las = 2
        )
        text(
            length(categories),
            0.9*max(dat[,5]),
            labels = sprintf(
                "Threshold :\n%3.2f",
                th_FoldChange
            ),
            col = "red"
        )
        boxplot(
            lfc,
            pch = 20,
            col = colors,
            main = "Fold Change of categories",
            ylab = "Fold Change [log2]",
            las = 2
        )
        abline(
            h = log2(th_FoldChange),
            col = "red",
            lty = "dashed"
        )

        print("Plot ranked genomic regions based on their Fold Change")
        plot(
            sort(
                dat[
                    -which(dat[,4]=="Random" | dat[,4]=="Negative"),
                    5
                ]
            ),
            pch = 20,
            main = "Activity of genomic regions\n(Random/Negative not included)",
            ylab = "Fold Change",
            xlab = "Ranked genomic regions"
        )
        idx <- which(dat[,5] >= th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
        abline(
            v = nrow(dat)-length(idx),
            col = "red",
            lty = "dashed"
        )
        text(
            nrow(dat)-length(idx),
            0.9*max(dat[,5]),
            labels = length(idx),
            pos = 2,
            offset = 0,
            col = "red"
        )
        idx <- which(
            dat[,5] < th_FoldChange &
            dat[,4] != "Random" &
            dat[,4] != "Negative"
        )
        text(
            length(idx)/2,
            0.9*max(dat[,5]),
            labels = length(idx),
            col = "black"
        )

        #- info sur l"analyse
        plot(
            c(0, 1),
            c(0, 1),
            ann = FALSE,
            bty = "n",
            type = "n",
            xaxt = "n",
            yaxt = "n"
        )
        text(
            x = 0.5,
            y = 1,
            "INFORMATION",
            cex = 1.5,
            pos = 1,
            offset = 0,
            col = "black"
        )
        text(
            x = 0,
            y = 0.70,
            labels = snakemake@wildcards[["id_sample"]],
            pos = 4,
            offset = 0,
            col = "black"
        )
        text(
            x = 0,
            y = 0.60,
            paste(
                "Method: FDR=",
                fdr,
                sep = ""
            ),
            pos = 4,
            offset = 0,
            col = "black"
        )
        dev.off()
        
        dat_groups <- data.frame(
            dat[,1:4],
            groups
        )
        write.table(
            dat_groups,
            file = paste0(
                snakemake@params[["outdir"]],
                "Groups.FDR=",
                fdr,
                ".grp"
            ),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t"
        )
    }
}


## Compute threshold base on inflexion point

#---  code from ROSE tool to dertermine super-enhancers  ---#
#--
#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
    inputVector <- sort(inputVector)
    inputVector[inputVector<0] <- 0 #set those regions with more control than ranking equal to zero
    slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
    xPt <- floor(
        optimize(
            numPts_below_line,
            lower = 1,
            upper = length(inputVector),
            myVector = inputVector,
            slope = slope
        )$minimum
    ) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
    y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.

    if(drawPlot){  #if TRUE, draw the plot
        plot(
            1:length(inputVector),
            inputVector,
            type = "l",
            ...
        )
        b <- y_cutoff-(slope* xPt)
        abline(
            v = xPt,
            h = y_cutoff,
            lty = 2,
            col = 8
        )
        points(
            xPt,
            y_cutoff,
            pch = 16,
            cex = 0.9,
            col = 2
        )
        abline(
            coef = c(b,slope),
            col = 2
        )
        title(
            paste(
                "x=",
                xPt,
                "\ny=",
                signif(
                    y_cutoff,
                    3
                ),
                "\nFold over Median=",
                signif(
                    y_cutoff/median(inputVector),
                    3
                ),
                "x\nFold over Mean=",
                signif(
                    y_cutoff/mean(inputVector),
                    3
                ),
                "x",
                sep=""
            )
        )
        axis(
            1,
            sum(inputVector==0),
            sum(inputVector==0),
            col.axis = "pink",
            col = "pink"
        ) #Number of regions with zero signal
    }
    return(
        list(
            absolute = y_cutoff,
            overMedian = y_cutoff/median(inputVector),
            overMean = y_cutoff/mean(inputVector)
        )
    )
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
    yPt <- myVector[x]
    b <- yPt-(slope*x)
    xPts <- 1:length(myVector)
    return(sum(myVector<=(xPts*slope+b)))
}
#--
#---  code from ROSE tool to determine super-enhancers  ---#

# determination du seuil
idx <- which(
    dat[,4] != "Negative" &
    dat[,4] != "Random"
)
th_FoldChange <- calculate_cutoff(dat[idx,5], drawPlot=F)$absolute

# identification des regions actives/inactives
idx <- which(dat[,5] >= th_FoldChange)
groups <- rep(
    "Inactive",
    nrow(dat)
)
groups[idx] <- "Active"
groups[which(is.na(dat[,5]))] <- "NA"

pdf(snakemake@output[["pdf"]])
par(mfrow=c(2,2))

#- boxplot en fonction des categories
categories <- unique(dat[,4])
fc <- list() ; lfc <- list()
for (category in categories) {
    idx <- which(dat[,4] == category)
    fc[[category]] <- dat[idx,5]
    lfc[[category]] <- log2(dat[idx,5])
}
colors <- rep("lightgrey",length(categories))
if ("Random" %in% categories) { colors[which(categories == "Random")] <- "darkblue" }
if ("PosEpromoter" %in% categories) { colors[which(categories == "PosEpromoter")] <- "darkgreen" }
if ("Positive" %in% categories) { colors[which(categories == "Positive")] <- "darkgreen" }
if ("Negative" %in% categories) { colors[which(categories == "Negative")] <- "darkred" }
boxplot(fc, pch=20, col=colors, main="Fold Change of categories", ylab="Fold Change", las=2)
text(length(categories), 0.9*max(dat[,5]), labels=sprintf("Threshold :\n%3.2f",th_FoldChange), col="red")
boxplot(lfc, pch=20, col=colors, main="Fold Change of categories", ylab="Fold Change [log2]", las=2)
abline(h=log2(th_FoldChange), col="red", lty="dashed")

#- ranked genomic regions based on their Fold Change
plot(sort(dat[which(dat[,4]!="Random" & dat[,4]!="Negative"),5]), pch=20, main="Activity of genomic regions\n(Random/Negative not included)", ylab="Fold Change", xlab="Ranked genomic regions")

idx <- which(dat[,5] >= th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
abline(v=nrow(dat)-length(idx), col="red", lty="dashed")
text(nrow(dat)-length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), pos=2, offset=0, col="red")
idx <- which(dat[,5] < th_FoldChange & dat[,4] != "Random" & dat[,4] != "Negative")
text(length(idx)/2, 0.9*max(dat[,5]), labels=length(idx), col="black")

print("debug_test4")

#- info sur l"analyse
plot(
    c(0, 1),
    c(0, 1),
    ann = FALSE,
    bty = "n",
    type = "n",
    xaxt = "n",
    yaxt = "n"
)
text(
    x = 0.5,
    y = 1,
    "INFORMATION",
    cex = 1.5,
    pos = 1,
    offset = 0,
    col = "black"
)
text(
    x = 0,
    y = 0.70,
    labels = snakemake@wildcards[["id_sample"]],
    pos = 4,
    offset = 0,
    col = "black"
)
text(
    x = 0,
    y = 0.60,
    "Method: Inflexion point",
    pos = 4,
    offset = 0,
    col = "black"
)
dev.off()

dat <- data.frame(
    dat[,1:4],
    groups
)
write.table(
    dat,
    file = snakemake@output[["groups"]],
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
)

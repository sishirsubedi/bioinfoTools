#library
library(circlize)
circos.clear()

df_genome <- read.csv("genomic_track.csv",header=1)
df_orf <- read.csv("orf_track.csv",header=1)
df_snp <- read.csv("snp_count_withtrack_plot.csv",header=1)

png("test.png",,width=2500,height=2500,res=300)
circos.genomicInitialize(df_orf)
circos.track(ylim = c(0, 0.5),
    bg.col = "black",
    bg.border = NA, track.height = 0.05)
dev.off()

circos.genomicInitialize(df_orf)
bgcol = rep(c("#5FFF4D", "#FF4D59"), 7)
circos.track(ylim = c(0, 1),
    bg.col = bgcol,
    bg.border = NA, track.height = 0.05)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        circos.lines(sort(runif(20)), runif(20), col = 4)
    })
circos.lines(x1, y1, pch = 16, cex = 0.5, type="h" , col="#69b3a2" , lwd=3)
dev.off()


df <- read.csv("snp_count_withtrack_plot.csv",header=1)
df = subset(df,df$sequence != '-')

png("test.png")


circos.initialize(factors = factors, xlim = c(0, 1))
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA )

# Now we can add a plot in this section! You can repeat these steps to add several regions
circos.lines(x1, y1, pch = 16, cex = 0.5, type="h" , col="#69b3a2" , lwd=3)

# Add axis
circos.axis(h="bottom" , labels.cex=0.4, direction = "inside" )

#clear
circos.clear()

dev.off()

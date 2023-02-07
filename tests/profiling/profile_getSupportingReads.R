library(profvis)
library(SVseek)

profvis({
    getSupportingReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
})

library(profvis)
library(SVseek)

profvis({
    getNonSupportingReads(data_fusion$gr, data_fusion$breakpoint, "qiaseq")
})

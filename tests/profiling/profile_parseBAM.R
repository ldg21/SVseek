library(profvis)
library(SVseek)

profvis({
    path <- system.file("extdata", package = "SVseek")
    input_bam <- file.path(path, "bams", "data_deletion_hg38.bam")
    parseBAM(input_bam)
})

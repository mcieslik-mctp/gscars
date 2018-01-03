

function() {
    countBarcodes("inst/extdata/test_2.fq", "B", "inst/extdata/4M-with-alts-february-2016.txt")

    countBarcodes("inst/extdata/test_2.fq", "B", "inst/extdata/4M-with-alts-february-2016.txt")    


    library(devtools)
    setwd("/mctp/users/mcieslik/proj/code/gscars/gscars")
    load_all()
    preprocessFastq("", "inst/extdata/test_1.fq", "inst/extdata/test_2.fq", "A", "B")
    
}

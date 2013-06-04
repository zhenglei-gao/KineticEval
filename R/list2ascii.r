##' Function to write a list into a file.
##'
##' @title Write the results
##' @param x a list.
##' @param file file name
##' @return NULL
##' @export
list2ascii <- function(x,file=paste(deparse(substitute(x)),".txt",sep="")) {
   # From the 'crosshybDetector' package from Paolo Uva <paolo.uva@gmail.com>
   # MHP July 7, 2004
   # R or S function to write an R list to an ASCII file.
   # This can be used to create files for those who want to use
   # a spreadsheet or other program on the data.
   #
   tmp.wid = getOption("width")  # save current width
   options(width=10000)          # increase output width
   sink(file)                    # redirect output to file
   print(x)                      # print the object
   sink()                        # cancel redirection
   options(width=tmp.wid)        # restore linewidth
   return(invisible(NULL))       # return (nothing) from function

}

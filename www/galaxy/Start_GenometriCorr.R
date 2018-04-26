# Start_GenometriCorr.R

###################################################
#                                                 #
#  command-line interface to GenometriCorr        #
#  functions, for use with Galaxy.                #
#                                                 #
###################################################

capture.output <- function (result, pdffile, output_options)
{
   if(output_options != "stats")
   {
      pdf(file=pdffile, width=10, height=19, paper="special")
   
      if (output_options != "vis")   #need to do a plot
      {
         mymat <- matrix(ncol=3, nrow=4)
         mymat[1,1] <- 1
         mymat[1,2] <- 2
         mymat[1,3] <- 3
         mymat[2,1] <- 4
         mymat[2,2] <- 5
         mymat[2,3] <- 6
         mymat[3,1] <- 7
         mymat[3,2] <- 8
         mymat[3,3] <- 9
         mymat[4,1] <- 10
         mymat[4,2] <- 11
         mymat[4,3] <- 12
       
         layout(mymat, heights=c(0.2,0.2,0.2,0.2))
         graphical.report(result, pdffile, make.new=FALSE)
      }
      if (output_options != "plot")  #need to do the bigger graphic
      {
         mymat <- matrix(ncol=2, nrow=8)
         mymat[1,1] <- 2
         mymat[1,2] <- 3
         mymat[2,1] <- 4
         mymat[2,2] <- 4
         mymat[3,1] <- 1
         mymat[3,2] <- 1
         mymat[4,1] <- 5
         mymat[4,2] <- 6
         mymat[5,1] <- 7
         mymat[5,2] <- 7
         mymat[6,1] <- 8
         mymat[6,2] <- 9 
         mymat[7,1] <- 10
         mymat[7,2] <- 10
         mymat[8,1] <- 11
         mymat[8,2] <- 12
         layoutresults <- 3
         
         layout(mymat, heights=c(0.05,0.05,0.15,0.15,0.15,0.15,0.15,0.15))
         visualize(result, pdffile, make.new=FALSE) 
      }
      dev.off()
   } 
   
   if (output_options == "stats")
   {
      show(result)
   }
}



# Reads the command line arguments
args <- commandArgs(trailingOnly=T)

suppressWarnings(suppressPackageStartupMessages(library('GenometriCorr',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('graphics',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('gdata',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('gplots',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('gtools',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('caTools',  warn.conflicts=F, verbose=F)))
suppressWarnings(suppressPackageStartupMessages(library('grid',  warn.conflicts=F, verbose=F)))



# Variables
query_file <- ""
reference_file <- ""
config_file <- ""
output_options <- ""

# Parse the command line arguments

config_file <- args[1]
query_file <- as.character(args[2])
reference_file <- as.character(args[3])
output_options <- args[4]
pdffile <- args[5]

conf<-new("GenometriCorrConfig",config_file)


result<-suppressWarnings(suppressPackageStartupMessages(GenometriCorr:::run.config(conf,query=query_file,reference=reference_file)))

hideoutput <- capture.output(result, pdffile=args[5], output_options)


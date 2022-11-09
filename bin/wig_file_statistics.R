filelist_fw <- c(
  "RVFV_virions/raw-wig-files/A_rep_1_forward.wig",
  "RVFV_virions/raw-wig-files/A_rep_2_forward.wig",
  "RVFV_virions/raw-wig-files/A_rep_3_forward.wig"
)

filelist_rev <- c(
  "RVFV_virions/raw-wig-files/A_rep_1_reverse.wig",
  "RVFV_virions/raw-wig-files/A_rep_2_reverse.wig",
  "RVFV_virions/raw-wig-files/A_rep_3_reverse.wig"
)

library(wig)

chromosome_lengths <- list(
  "",
  "",
  ""
)

import_wig_files<-function(files=vector())
{
  chromosomes = vector();
  chromosome_length = list();
  chromosome_values = list();
  chr_dat = list();
  filelist = list();
  
  for(i in seq(1,length(files)))
  {
    filelist[[i]] <- import_wig(files[i])
    # generate a list of chromosome names
    chromosomes=levels(factor(c(chromosomes, levels(as.factor(filelist[[i]]$chr)))))
  }
  
  # estimate the length of the chromosomes
  for(chr in chromosomes)
  {
    chromosome_length[[chr]] = 0
    for (i in seq(1,length(files)))
    {
      chromosome_length[[chr]] = max(
        chromosome_length[[chr]],
        as.vector(filelist[[i]][filelist[[i]]$chr==chr,]$pos)
      );
    }
  }
  
  # generate a list of chromosome matrices with all values
  for(chr in chromosomes)
  {
    chromosome_values[[chr]] = matrix(nrow=chromosome_length[[chr]], ncol=length(files),data=c(0))
    # loop through our input files and import the values into our matrix
    for (i in seq(1,length(files)))
    {
      lines     = filelist[[i]]$chr==chr
      positions = filelist[[i]][lines,]$pos
      values    = filelist[[i]][lines,]$val
      chromosome_values[[chr]][positions, i] = values
    }
  }
  
  # build the result list
  for(chr in chromosomes)
  {
    chr_dat[[chr]] = list(name=chr, len=chromosome_length[[chr]], val=chromosome_values[[chr]])
  }
  
  return(chr_dat);
}

do_cor_analysis<-function(files=vector(), dat=list(), output_filename_dir=paste(getwd(),"/",sep=""))
{
  for(chr in seq(1, length(dat)))
  {
    chr_name <- dat[[chr]]$name
    chr_len <- dat[[chr]]$len
    cat("Chromosome: ", chr_name, "; length ", chr_len, " bp\n");
    
    cor_res <- cor(dat[[chr]]$val)
    colnames(cor_res)<-basename(files)
    rownames(cor_res)<-colnames(cor_res)
    
    csv_output_file<-paste(output_filename_dir, chr_name, "_", chr_len, "bp.csv", sep="")
    cat("Writing tabular output to ", csv_output_file, "...\n");
    write.csv(cor_res, file=csv_output_file)
    
    print(cor_res)
    
    # print all the stuff
    heatmap(cor_res,symm=TRUE)
    heatmap_output_file<-paste(output_filename_dir, chr_name, "_", chr_len, "bp.pdf", sep="")
    pdf(file=heatmap_output_file, paper="a4")
    heatmap(cor_res,symm=TRUE)
    dev.off()
  }
}

do_cor_analysis(filelist_fw, import_wig_files(filelist_fw), output_filename_dir = paste(getwd(), "/forward_", sep=""))
do_cor_analysis(filelist_rev, import_wig_files(filelist_rev), output_filename_dir = paste(getwd(), "/reverse_", sep=""))
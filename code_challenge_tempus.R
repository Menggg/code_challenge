#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)#first argument is vcf file name, second argyment is gtf/gff file name
print(paste("argument1 is", args[1]))
print(paste("argument2 is", args[2]))
##loading library, functions and data sets
#install.packages('stringr')
library(stringr)
larsub<-function(x) {
  #input is two string vector
  #output is overlaped strings vector between two input strings
  str_a<-x[1]
  str_b<-x[2]
  prior=NULL
  for(j in 1:nchar(str_a)){
    sb<-unique(combn(strsplit(str_a, "")[[1]],j, FUN=paste, collapse=""))
    if(length(unlist(str_extract_all(str_b,sb)))==0){ 
      post<-prior
      return(post)
    }
    prior<-unlist(str_extract_all(str_b,sb))
  }
}
sub_chr <- function(x){
  x=gsub("chr","",x,fixed=TRUE)
}
sub_af <- function(x){
  x=gsub("allele_freq: ","",x,fixed=TRUE)
}
sub_gene <- function(x){
  x=gsub("Gene: ","",x,fixed=TRUE)
}
sub_var <- function(x){
  x=gsub("Existing_variation: ","",x,fixed=TRUE)
}


####reading input data
input=read.table(args[1],stringsAsFactors=FALSE)
input=as.data.frame(input[,c(1:9,3,10,3,11,3,3)])
print("===VCF file is loaded ===")
#reading gff file and get gene ID
print("===reading gtf file ===")
gtf=read.delim(args[2], head=F, comment.char="#")
print("===gtf file is loaded ===")
gff.genes <- gtf[,1:5]
chr=gtf[,1]
chr<- sub_chr(chr)
gff.genes[,1]=chr
base="http://exac.hms.harvard.edu/"
item="/rest/variant/"

##start loop all the variants
print("===looping variants===")
for(n in 1:nrow(input)){
  ####pasing and identify type of mutation
  term=''
  temp=''
  mut_length_r = str_length(input[n,4])
  mut_length_a = str_length(input[n,5])
  n_mut_r = str_split(input[n,4],'\\,')[[1]]
  n_mut_a = str_split(input[n,5],'\\,')[[1]]
  if(mut_length_r == 1 & mut_length_a == 1){
    #type is substitution
    term = paste(term,"substitution",sep='')
    pos_mut = input[n,2]
    #locate mutation in gff
    loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut & as.numeric(gff.genes[,5]) >= pos_mut)
    loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
    loc=loc_1[loc_2]
    temp_gff=gff.genes[loc,]
    n_exon = length(grep("exon",temp_gff[,3]))
    n_gene = length(grep("gene",temp_gff[,3]))
    n_all=length(temp_gff[,3])
    if(n_all<1){
      #back intergenic
      term = paste(term,";intergenic",sep='')
    }else if(n_exon>0){
      #back exon
      term = paste(term,";exon",sep='')
    }else if(n_gene>0){
      #back intron
      term = paste(term,";intron",sep='')
    }
  }else if(length(n_mut_a)==1){
    if(mut_length_r<10 & mut_length_a <10){
      #can't deal with long segment due to computing speed 
      overlap1 = larsub(c(as.character(input[n,4]),as.character(input[n,5])))
      overlap2 = larsub(c(as.character(input[n,5]),as.character(input[n,4])))
    }
    
    diff_length = abs(mut_length_r-mut_length_a)
    if(mut_length_r<mut_length_a){
      #type is insertion
      term = paste(term,"insertion",sep='')
      if(mut_length_r<10 & mut_length_a <10){
        mut_offset=0
        if(length(overlap1)>0)
          mut_offset = str_locate(as.character(input[n,5]),overlap1[1])[2]
        pos_mut = mut_offset+input[n,2]
        #locate mutation in gff
        loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut & as.numeric(gff.genes[,5]) >= pos_mut)
        loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
        loc=loc_1[loc_2]
        temp_gff=gff.genes[loc,]
        n_exon = length(grep("exon",temp_gff[,3]))
        n_gene = length(grep("gene",temp_gff[,3]))
        n_all=length(temp_gff[,3])
        if(n_all<1){
          #back intergenic
          term = paste(term,";intergenic",sep='')
        }else if(n_exon>0){
          #back exon
          term = paste(term,";exon",sep='')
        }else if(n_gene>0){
          #back intron
          term = paste(term,";intron",sep='')
        }
      }
    }
    
    if(mut_length_r>mut_length_a){
      #type is deletion
      term = paste(term,"deletion",sep='')
      if(mut_length_r<10 & mut_length_a <10){
        mut_offset=0
        if(length(overlap2)>0)
        mut_offset = str_locate(as.character(input[n,4]),overlap2[1])[2]+1
        pos_mut = mut_offset+input[n,2]
        #locate mutation in gff
        loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut & as.numeric(gff.genes[,5]) >= pos_mut)
        loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
        loc=loc_1[loc_2]
        temp_gff=gff.genes[loc,]
        n_exon = length(grep("exon",temp_gff[,3]))
        n_gene = length(grep("gene",temp_gff[,3]))
        n_all=length(temp_gff[,3])
        if(n_all<1){
          #back intergenic
          term = paste(term,";intergenic",sep='')
        }else if(n_exon>0){
          #back exon
          term = paste(term,";exon",sep='')
        }else if(n_gene>0){
          #back intron
          term = paste(term,";intron",sep='')
        }
      }
    }
    if(mut_length_r==mut_length_a){
      #type is substitution
      term = paste(term,"substitution",sep='')
      if(mut_length_r<10 & mut_length_a <10){
        if(length(overlap1)>0)
          mut_offset = str_locate(as.character(input[n,4]),overlap1[1])
        if(length(overlap1)==0)
          mut_offset = c(2,str_length(input[n,4]))
        pos_mut = mut_offset[1]-1+input[n,2]
        pos_mut = c(pos_mut,mut_offset[2]+input[n,2])
        pos_mut = pos_mut[which(pos_mut>input[n,2])]
        
        for(i in 1:length(pos_mut)){
          #locate mutation in gff
          loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut[i] & as.numeric(gff.genes[,5]) >= pos_mut[i])
          loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
          loc=loc_1[loc_2]
          temp_gff=gff.genes[loc,]
          n_exon = length(grep("exon",temp_gff[,3]))
          n_gene = length(grep("gene",temp_gff[,3]))
          n_all=length(temp_gff[,3])
          if(n_all<1){
            #back intergenic
            term = paste(term,";intergenic",sep='')
          }else if(n_exon>0){
            #back exon
            term = paste(term,";exon",sep='')
            temp = term
          }else if(n_gene>0){
            #back intron
            term = paste(term,";intron",sep='')
          }
          term=temp
      }
      }
    }
  }else if(length(n_mut_a)>1){
    for(i in 1:length(n_mut_a)){
      term=''
      temp=''
      if(mut_length_r<10 & mut_length_a <10){
        #can't deal with long segment due to computing speed
        overlap1 = larsub(c(as.character(input[n,4]),as.character(n_mut_a[i])))
        overlap2 = larsub(c(as.character(n_mut_a[i]),as.character(input[n,4])))
        mut_offset = abs(mut_length_r-mut_length_a)
      }
      
      if(mut_length_r<mut_length_a){
        #type is insertion
        term = paste(term,"insertion",sep='')
        if(mut_length_r<10 & mut_length_a <10){
          mut_offset=0
          if(length(overlap1)>0)
          mut_offset = str_locate(as.character(n_mut_a[i]),overlap1[1])[2]
          pos_mut = mut_offset+input[n,2]
          #locate mutation in gff
          loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut & as.numeric(gff.genes[,5]) >= pos_mut)
          loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
          loc=loc_1[loc_2]
          temp_gff=gff.genes[loc,]
          n_exon = length(grep("exon",temp_gff[,3]))
          n_gene = length(grep("gene",temp_gff[,3]))
          n_all=length(temp_gff[,3])
          if(n_all<1){
            #back intergenic
            term = paste(term,";intergenic",sep='')
          }else if(n_exon>0){
            #back exon
            term = paste(term,";exon",sep='')
            temp = term
          }else if(n_gene>0){
            #back intron
            term = paste(term,";intron",sep='')
          }
        }
      }
      
      if(mut_length_r>mut_length_a){
        #type is deletion
        term = paste(term,"deletion",sep='')
        if(mut_length_r<10 & mut_length_a <10){
          mut_offset=0
          if(length(overlap1)>0 & length(overlap2)>0)
          mut_offset = str_locate(as.character(input[n,4]),overlap2[1])[2]+1
          pos_mut = mut_offset+input[n,2]
          #locate mutation in gff
          loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut & as.numeric(gff.genes[,5]) >= pos_mut)
          loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
          loc=loc_1[loc_2]
          temp_gff=gff.genes[loc,]
          n_exon = length(grep("exon",temp_gff[,3]))
          n_gene = length(grep("gene",temp_gff[,3]))
          n_all=length(temp_gff[,3])
          if(n_all<1){
            #back intergenic
            term = paste(term,";intergenic",sep='')
          }else if(n_exon>0){
            #back exon
            term = paste(term,";exon",sep='')
            temp = term
          }else if(n_gene>0){
            #back intron
            term = paste(term,";intron",sep='')
          }
        }
      }
      if(mut_length_r==mut_length_a){
        #type is substitution
        term = paste(term,"substitution",sep='')
        if(mut_length_r<10 & mut_length_a <10){
          if(length(overlap1)>0)
            mut_offset = str_locate(as.character(input[n,4]),overlap1[1])
          if(length(overlap1)==0)
            mut_offset = c(2,str_length(input[n,4]))
          pos_mut = mut_offset[1]-1+input[n,2]
          pos_mut = c(pos_mut,mut_offset[2]+input[n,2])
          pos_mut = pos_mut[which(pos_mut>input[n,2])]
          
          for(i in 1:length(pos_mut)){
            #locate mutation in gff
            loc_1=which(as.numeric(gff.genes[,4]) <= pos_mut[i] & as.numeric(gff.genes[,5]) >= pos_mut[i])
            loc_2=which(as.character(gff.genes[loc_1,1]) == as.character(input[n,1]))
            loc=loc_1[loc_2]
            temp_gff=gff.genes[loc,]
            n_exon = length(grep("exon",temp_gff[,3]))
            n_gene = length(grep("gene",temp_gff[,3]))
            n_all=length(temp_gff[,3])
            if(n_all<1){
              #back intergenic
              term = paste(term,";intergenic",sep='')
            }else if(n_exon>0){
              #back exon
              term = paste(term,";exon",sep='')
              temp = term
            }else if(n_gene>0){
              #back intron
              term = paste(term,";intron",sep='')
            }
          }
          term=temp
        }
        
      }
    }
  }
  
  if(str_length(term)>0)
    input[n,10]=term
  ##extract DP, AO and AO/(RO+AO)
  term=''
  reads_count = str_split(input[n,11],"\\:")[[1]]
  DP = reads_count[3]
  AO = reads_count[5]
  if(length(n_mut_a)>1){
    reads_r = str_split(reads_count[7],"\\,")[[1]]
    r_AO1 = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_r[1])),3)
    r_AO2 = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_r[2])),3)
    r_AO = paste(r_AO1,r_AO2,sep=',')
  }else{
    r_AO = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_count[7])),3)
  }
  term = paste(DP,AO,r_AO,sep=';')
  if(str_length(term)>0)
    input[n,12]=term
  term=''
  reads_count = str_split(input[n,13],"\\:")[[1]]
  DP = reads_count[3]
  AO = reads_count[5]
  if(length(n_mut_a)>1){
    reads_r = str_split(reads_count[7],"\\,")[[1]]
    r_AO1 = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_r[1])),3)
    r_AO2 = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_r[2])),3)
    r_AO = paste(r_AO1,r_AO2,sep=',')
  }else{
    r_AO = round(as.numeric(reads_count[5])/(as.numeric(reads_count[5])+as.numeric(reads_count[7])),3)
  }
  term = paste(DP,AO,r_AO,sep=';')
  if(str_length(term)>0)
    input[n,14]=term
  ##extract info from ExAC
  term = ''
  variant=paste(input[n,1],input[n,2],input[n,4],input[n,5],sep='-')
  URL_term=paste(base,item,variant,sep='')
  system(paste("curl",URL_term,"> json"),ignore.stderr = T)
  text_json = suppressWarnings(read.delim("json",head=F))
  jfile = as.data.frame(str_split(text_json[1,1],"\\, "))
  af=as.character(jfile[grep("allele_freq:",jfile[,1])[1],1])
  af=sub_af(af)
  gene_ID=as.character(jfile[grep("Gene:",jfile[,1])[1],1])
  gene_ID=sub_gene(gene_ID)
  var_ID=as.character(jfile[grep("Existing_variation:",jfile[,1])[1],1])
  var_ID=sub_var(var_ID)
  term = paste(af,gene_ID,var_ID,sep=';')
  if(str_length(term)>0)
    input[n,15]=term

  print(n)
}
##export to local folder
colnames(input) = c("#CHROM","POS",	"ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TYPE","normal","DP/AO/r_AO","vaf5","DP/AO/r_AO","info_ExAC")
write.table(input,"code_challege_annoted.vcf",col.names = T,row.names = F,quote=F,na="NA")
print("====It is done====")

require(plyr)
require(Biostrings)

## program...

extractAcc = function (acc, upstream, downstream) {
  if (grepl("[A-Z]{3}[0-9]{5}\\.[0-9]{1}",acc)) {
    prot = F
  } else {
    prot = T
    if (grepl("_",acc)) {
      acc <- sub("(.*)_.*","\\1",acc)
    }
    url = paste("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=",acc,"&format=default&style=raw&Retrieve=Retrieve",sep="")
    page = paste(readLines(url),collapse="")
    if (grepl("ERROR 12 No entries found.",page)) {
      return("obsolete")
      } else {
      #get metadata on strain, NCBI taxid, ENA accession number and the complete aa translation through regex
      strain = sub(".+;(.*?)OC\\s.+", "\\1", page)
      strain = gsub("OS\\s{2}", "", strain)
      strain = sub("\\s", "", strain)
      if (grepl("OG",strain)) {
        OG = strsplit(strain,"OG")
        strain = OG[[1]][1]
        organelle = OG[[1]][2:length(OG[[1]])]
        organelle = sub("\\s{1,}", "", organelle)
        organelle = paste(organelle,collapse=" ")
      } else {
        organelle = "genomic"
        }
      if (grepl("DOI=",page)) {
        doi = sub(".+DOI=(.*?);.+", "\\1", page)
      } else {
          doi = "NA"
          }
      NCBI.id = sub(".+NCBI_TaxID=([0-9]*).+", "\\1", page)
      #if (grepl(";",NCBI.id)) {NCBI.id = strsplit(NCBI.id, ";")[[1]][1]}
      EMBL.id = sub(".+EMBL;.*([A-Z]{3}[0-9]{5}\\.[0-9]{1});.+", "\\1", page)
      #if (nchar(EMBL.id)>12) {EMBL.id = sub(".+EMBL;.*?\\s{1,}(.*?\\.[0-9]);.+", "\\1", page)}
      aa.seq = sub(".+CRC64;(.*?)//", "\\1", page)
      #remove whitespace from the aa sequence
      aa.seq = gsub("[[:space:]]", "", aa.seq)
      #return data frame containing all the extracted information
      outtable = data.frame(strain=strain,origin=organelle,DOI=doi,NCBI.id=NCBI.id,EMBL.id=EMBL.id,UniProtKB.id=acc,aa.seq=aa.seq)
      acc = outtable$EMBL.id
    }
  }
    #read in the page record from ena
    url <- paste("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence&id=",acc,"&format=default&style=raw&Retrieve=Retrieve",sep="")
    page <- readLines(url)
    
    #get the alternative accession number at this point for the genomic DNA record
    alt_id <- grep("ID\\s{3}",page,value=T)
    alt_id <- sub("ID\\s{3}(.*?);.+","\\1",alt_id)
    #take the UniProtKb/TrEMBL id value
    id <- grep("TrEMBL:",page,value=T)
    #in rare cases the record has more than one protein entry, the following statement sorts this out
    if (length(id)>1) {
      #find which entry contains our wanted protein and extract only that region from the record
      match = grep(acc,page)
      CDS_match = grep("FT   CDS",page)
      dis = match-CDS_match
      region_start = match-min(dis[which(dis>0)])
      page = page[region_start:match]
      id = grep("TrEMBL:",page,value=T)
    }
    id <- sub(".+TrEMBL:(.*?)[[:punct:]]","\\1",id)
    #if the accession number is not directly supplied here, get it from uniprot
    if (length(id)==0) {
      url_uniprot = paste("http://www.uniprot.org/uniprot/?query=",acc,"&columns=id&format=tab",sep="")
      uniprot <- readLines(url_uniprot)
      id <- uniprot[2]
      #or if not found, from uniparc?
      if (is.na(id)) {
        acc_uniparc <- strsplit(as.character(acc),"\\.")[[1]][1]
        url_uniparc <- paste("http://www.uniprot.org/uniparc/?query=",acc_uniparc,"&format=tab",sep="")
        uniparc <- read.table(url_uniparc,sep="\t",header=T)
        id_split <- strsplit(as.character(uniparc$UniProtKB),"; ")
        #and grab the first non-obsolete record for "id"
        id <- id_split[[1]][!grepl("obsolete",id_split[[1]])][1]
        #if this doesn't return a valid accession number, mark the entry for removal as redundant
        if (is.na(id)) {
          id <- "obsolete"
        }
      }
    }
    if (id=="obsolete") {
      return("obsolete")
      } else {
    #get information on the genomic region to extract from (nucleotides)
    #and expand the range according to the input to the function
    CDS <- grep("FT   CDS",page,value=T)
    #sometimes the record is joined, so check for that and run the alternative pipeline (for comments, check the "else" part)
    if (grepl("join",CDS)) {
      CDS = gsub("<|>","",CDS)
      CDS2 <- page[grep("FT   CDS",page)+1]
      if (grepl("\\.{2}",CDS2)==F) {
        split_CDS = strsplit(CDS,",")
        CDS = split_CDS[[1]][1]
        CDS2 = split_CDS[[1]][2]
      } else {
        CDS = strsplit(CDS,",")[[1]]
      }
      CDS2 = gsub("<|>","",CDS2)
      if (grepl(":",CDS)) {
        from = as.numeric(sub(".+:([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS))
      } else {
        from = as.numeric(sub(".+\\(([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS))
      }
      if (grepl(":",CDS2)) {
        from_CDS2 = as.numeric(sub(".+:([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS2))
      } else {
        from_CDS2 = as.numeric(sub(".+\\(([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS2))
      }
      genome = sub(".+join\\((.*?):.+", "\\1",CDS)
      to = as.numeric(sub(".+\\.{2}([0-9]*?)", "\\1", CDS))
      to_CDS2 = as.numeric(sub(".+\\.{2}([0-9]*?)\\)+", "\\1", CDS2))
      to2 = to
      from2_CDS2 = from_CDS2
      if (grepl("complement",CDS)) {
        from2 = from-downstream
        to2_CDS2 = to_CDS2+upstream
      } else {
        from2 = from-upstream
        to2_CDS2 = to_CDS2+downstream
      }
      if (grepl(":",CDS)==F) {
        genome = alt_id
        }
      url2 = paste("http://www.ebi.ac.uk/ena/data/view/",genome,"&display=fasta&range=",from2,"-",to2,sep="")
      out = readLines(url2, warn=F)
      out[2] = paste(out[2:length(out)],collapse="")
      out = out[1:2]
      url2_CDS2 = paste("http://www.ebi.ac.uk/ena/data/view/",genome,"&display=fasta&range=",from2_CDS2,"-",to2_CDS2,sep="")
      out_CDS2 = readLines(url2_CDS2, warn=F)
      if (grepl("complement",CDS)) {
        out[2] = as.character(reverseComplement(DNAString(x=out[2], start=1, nchar=NA)))
        out_CDS2 = paste(out_CDS2[2:length(out_CDS2)],collapse="")
        out_CDS2 = as.character(reverseComplement(DNAString(x=out_CDS2, start=1, nchar=NA)))
        out[2] = paste(out_CDS2,out[2],collapse="",sep="")
      } else {
        out[2] = paste(paste(out_CDS2[2:length(out_CDS2)],collapse=""),out[2],collapse="",sep="")
      }
      len = (to-from+1)+(to_CDS2-from_CDS2+1)
      if (nchar(out[2])>len) {
        out[2] = substr(out[2],from2,to2)
        }
      if (grepl("complement",CDS)) {
        out[1] = sub("Location.+",paste("Protein",acc,"CDS (complement)"),out[1])
        if (upstream!=0) {
          if (nchar(out[2])==len+upstream+downstream) {
            out[1] = paste(out[1],"-",upstream,"nt")
          } else {
              out[1] = paste(out[1],"-",upstream-((len+upstream+downstream)-nchar(out[2])),"nt")
          }
          }
        if (downstream!=0) {
          if (from2>1) {
            out[1] = paste(out[1],"+",downstream,"nt")
          } else {
              out[1] = paste(out[1],"+",from-from2,"nt")
          }
          }
      } else {
        out[1] = sub("Location.+",paste("Protein",acc,"CDS"),out[1]) 
      if (upstream!=0) {
        if (from2>1) {
          out[1] = paste(out[1],"-",upstream,"nt")
        } else {
            out[1] = paste(out[1],"-",from-from2,"nt")
        }
        }
      if (downstream!=0) {
        if (from2>1) {
          out[1] = paste(out[1],"+",(nchar(out[2])-len)-(upstream),"nt")
        } else {
            out[1] = paste(out[1],"+",(nchar(out[2])-len)-(from-from2),"nt")
        }
        }
      }
    } else {
      #sometimes "<" and ">" are present in the CDS line, so the next line removes them and saves us from a lot of trouble
      CDS = gsub("<|>","",CDS)
      #sometimes the CDS is in a different format - check and correct for this if necessary
      if (grepl(":",CDS)) {
        from = as.numeric(sub(".+:([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS))
      } else {
        from = as.numeric(sub(".+\\(([0-9]*?)[[:punct:]]{2}.+", "\\1", CDS))
      }
      #get the proper from and to, depending on if the region is in the complement or not
      if (grepl("complement",CDS)) {
        genome = sub(".+complement\\((.*?):.+", "\\1",CDS)
        to = as.numeric(sub(".+[[:punct:]]{2}([0-9]*?)\\)", "\\1", CDS))
        from2 = from-downstream
        to2 = to+upstream
      } else {
        genome = sub(".+CDS.*\\s(.*?):.+", "\\1",CDS)
        to = as.numeric(sub(".+\\.{2}([0-9]*?)", "\\1", CDS))
        from2 = from-upstream
        to2 = to+downstream
      }
      #if negative values of "from" appear, they're not retrieved properly from the database, so in this case the value is set to 1
      if (from2 < 1) {
        from2 = 1
        }
      #check and correct for proper accession number
      if (grepl(":",CDS)==F) {
        genome = alt_id
        }
      #the sequence itself is extracted in fasta format from the ENA database
      url2 = paste("http://www.ebi.ac.uk/ena/data/view/",genome,"&display=fasta&range=",from2,"-",to2,sep="")
      out = readLines(url2, warn=F)
      out[2] = paste(out[2:length(out)],collapse="")
      out = out[1:2]
      len = to-from+1
      #check if the record was taken out as requested
      if (nchar(out[2])>len) {
        out[2] = substr(out[2],from2,to2)
        }
      #sometimes seems that the first character of the entry is missed -> take out a longer region and remove the added characters
      if (nchar(out[2])==len-1) {
        url2 = paste("http://www.ebi.ac.uk/ena/data/view/",genome,"&display=fasta&range=",from2-1,"-",to2,sep="")
        out = readLines(url2, warn=F)
        out[2] = paste(out[2:length(out)],collapse="")
        out = out[1:2]
        out[2] = substr(out[2],2,nchar(out[2]))
      }
      #some R magic to obtain the sequence in the correct orientation including the defined upstream and downstream parts
      if (grepl("complement",CDS)) {
        out[1] = sub("Location.+",paste("Protein",acc,"CDS (complement)"),out[1])
        out[2] = as.character(reverseComplement(DNAString(x=out[2], start=1, nchar=NA)))
        if (upstream!=0) {
          if (nchar(out[2])==len+upstream+downstream) {
            out[1] = paste(out[1],"-",upstream,"nt")
          } else {
              out[1] = paste(out[1],"-",upstream-((len+upstream+downstream)-nchar(out[2])),"nt")
          }
          }
        if (downstream!=0) {
          if (from2>1) {
            out[1] = paste(out[1],"+",downstream,"nt")
          } else {
              out[1] = paste(out[1],"+",from-from2,"nt")
          }
          }
      } else {
        out[1] = sub("Location.+",paste("Protein",acc,"CDS"),out[1]) 
      if (upstream!=0) {
        if (from2>1) {
          out[1] = paste(out[1],"-",upstream,"nt")
        } else {
            out[1] = paste(out[1],"-",from-from2,"nt")
        }
        }
      if (downstream!=0) {
        if (from2>1) {
          out[1] = paste(out[1],"+",(nchar(out[2])-len)-(upstream),"nt")
        } else {
            out[1] = paste(out[1],"+",(nchar(out[2])-len)-(from-from2),"nt")
        }
        }
      }
      #^-- works perfectly for at least non-complementary sequences, CHECK FUNCTION FOR COMPLEMENTARY SHORT SEQUENCES (when going out of range)!
    }
    
    
    #then the table
    #make the correct url address by grabbing the UniProtKB/TrEMBL "id" from the beginning of the function
    if (prot==F) {
      url = paste("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=",id,"&format=default&style=raw&Retrieve=Retrieve",sep="")
      page = paste(readLines(url),collapse="")
      if (grepl("ERROR 12 No entries found.",page)) {
        url_uniprot = paste("http://www.uniprot.org/uniprot/?query=",acc,"&columns=id&format=tab",sep="")
        uniprot <- readLines(url_uniprot)
        id <- uniprot[2]
        #or if not found, from uniparc?
        if (is.na(id)) {
          acc_uniparc <- strsplit(as.character(acc),"\\.")[[1]][1]
          url_uniparc <- paste("http://www.uniprot.org/uniparc/?query=",acc_uniparc,"&format=tab",sep="")
          uniparc <- read.table(url_uniparc,sep="\t",header=T)
          id_split <- strsplit(as.character(uniparc$UniProtKB),"; ")
          #and grab the first non-obsolete record for "id"
          id <- id_split[[1]][!grepl("obsolete",id_split[[1]])][1]
          #if this doesn't return a valid accession number, mark the entry for removal as redundant
          if (is.na(id)) {
            id <- "obsolete"
          }
        }
        url = paste("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=",id,"&format=default&style=raw&Retrieve=Retrieve",sep="")
        page = paste(readLines(url),collapse="")
      }
      #get metadata on strain, NCBI taxid, ENA accession number and the complete aa translation through regex
      strain = sub(".+;(.*?)OC\\s.+", "\\1", page)
      strain = gsub("OS\\s{2}", "", strain)
      strain = sub("\\s", "", strain)
      if (grepl("OG",strain)) {
        OG = strsplit(strain,"OG")
        strain = OG[[1]][1]
        organelle = OG[[1]][2:length(OG[[1]])]
        organelle = sub("\\s{1,}", "", organelle)
        organelle = paste(organelle,collapse=" ")
      } else {
        organelle = "genomic"
        }
      if (grepl("DOI=",page)) {
        doi = sub(".+DOI=(.*?);.+", "\\1", page)
      } else {
          doi = "NA"
          }
      NCBI.id = sub(".+NCBI_TaxID=([0-9]*).+", "\\1", page)
      aa.seq = sub(".+CRC64;(.*?)//", "\\1", page)
      #remove whitespace from the aa sequence
      aa.seq = gsub("[[:space:]]", "", aa.seq)
      #return data frame containing all the extracted information
      outtable = data.frame(strain=strain,origin=organelle,DOI=doi,NCBI.id=NCBI.id,EMBL.id=acc,UniProtKB.id=id,aa.seq=aa.seq)
    }
    
    #format the headers to make further processing easier
    header_split = strsplit(out[1],"\\|")
    out[1] = paste(header_split[[1]][1],id,header_split[[1]][3],sep="|")
    out[1] = gsub("[[:space:]]","_",out[1])
    out[1] = gsub(",","",out[1])
    
    return <- list(out,outtable)
    return(return)
      }
    }


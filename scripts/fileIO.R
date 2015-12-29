#########################################################################################
# These are library functions related to IO of BioNano data files, such as map/xmap etc.
#
# Last modified 1/17/2014
# Direct questions/feedback to Heng Dai 
# hdai@bionanogenomics.com
##########################################################################################
 
read.bionano.header = function(f=NULL) {
# This function will read in the header section of a generic bioanano data file, like xmap etc
# The header is defined as leading lines with start with # or complete blank line as a continuation.
#   rule exception: to be back compatible with earlier version of .map file, header can also be "Software version: $.+$"
# The result is organized in a simple list with $headerLineCount and $headerLineArray

    #f = "W:/Research/Computational Biology/Dai/Arabidopsis_Merge/GenomeMap2-to-colv2_2.xmap"
    
    # we initiate result object
    # headerLineArray contains lines of header sections and headerLineCount is simply the count of this array
    # headerVer is defined such as XMAP Version:\t1.0 
    # headerCols is the array contains the col names as defined in #h, delimited by white space
    # headerColTypes is the class type of columns using int,float,string
    result = list(error = NULL, errorMsg = "", headerLineCount = 0, headerLineArray=NULL, headerVer = NULL, headerCols = NULL, headerColTypes = NULL)
    
    # we first check for file existence and access rights etc.
    if(!file.exists(f)) { 
        result$error = -100
        result$errorMsg = paste("Unable to read in file '", f , "'", sep="")
        return(result)
    } 
    
    # Now we open file handler
    fcon = file(f, 'r') 
    
    # we read one line at a time.
    lineIsHeader = TRUE
    headerLines = c()
    while (lineIsHeader) {
        nextLine = readLines(fcon, n=1)
        nextLine = sub("^\\s+", "", nextLine) # remove Leading white space
        
        if (substring(nextLine,1,1) == "#" | nchar(nextLine) == 0) {
            lineIsHeader = TRUE
        } else {
            # fix 
            if (length(grep("^Software version:",  nextLine, ignore.case = TRUE, perl=TRUE)) > 0) {
                lineIsHeader = TRUE
            } else {
                lineIsHeader = FALSE
            }
        }
        
        if (lineIsHeader) {
            headerLines = c(headerLines, nextLine)
            if (substring(nextLine,1,2) == "#h") { #h
                tmpS =  substring(nextLine,3,nchar(nextLine))
                tmpS = sub("^\\s+", "", tmpS)
                result$headerCols = strsplit(tmpS, "\\s+", perl=TRUE)[[1]]
            }
            
            if (substring(nextLine,1,2) == "#f") { #f
                tmpS =  substring(nextLine,3,nchar(nextLine))
                tmpS = sub("^\\s+", "", tmpS)
                result$headerColTypes = as.factor(strsplit(tmpS, "\\s+", perl=TRUE)[[1]])
            }
            
            if (length(grep("File Version:x\\s+.+$", nextLine, ignore.case = TRUE, perl=TRUE)) > 0) {
                result$headerVer = strsplit(nextLine, ":\\s+", perl=TRUE)[[1]][2]
            }
        } else {
            lineIsHeader = FALSE
            close(fcon)
        }
    }
    
    result$headerLineArray = headerLines
    result$headerLineCount = length(headerLines)
    
    return(result)
    
}

read.bionano.map = function(f=NULL) {
# This is the engine to read in the xmap/cmap/smap etc file as a data frame, usually it is called by read.xmap or read.cmap

    #f = "W:/Research/Computational Biology/Dai/Arabidopsis_Merge/GenomeMap2-to-colv2_2.xmap"

    # we initiate result object
    result = list(error = NULL, errorMsg = "", xmap = NULL)
   
    # first we read in header lines
    header.obj = read.bionano.header(f)
    if (!is.null(header.obj$error)) { 
        result$error = -100
        result$errorMsg = paste("Unable to read in file '", f , "'", sep="")
        return(result)
    } 
    
    # we will read in xmap as data frame
	m = as.data.frame(read.table(f, header=F, skip = header.obj$headerLineCount, fill=T))
    if (length(header.obj$headerCols) >= length(colnames(m))) {
        colnames(m) = header.obj$headerCols[1:NCOL(m)]
    } else {
        colnames(m)[1:NCOL(m)] = header.obj$headerCols
    }
    
    #TBD, class assignment
    
    return(m)
    # Test
    # read.bionano.map(f)
}

read.xmap = function(f=NULL, return.header=FALSE) {
# This function will read in the xmap file as a xmap data frame

    #f = "W:/Research/Computational Biology/Dai/Arabidopsis_Merge/GenomeMap2-to-colv2_2.xmap"
    m = read.bionano.map(f=f)
    # bug fix in case some time col name is RefcontigID in some earlier verisons
    if ("RefcontigID" %in% colnames(m)) {
        colnames(m)[colnames(m) == "RefcontigID"] = "RefContigID"
    }
    return(m)
}

read.cmap = function(f=NULL) {
# This function will read in the cmap file as a cmap data frame

    m = read.bionano.map(f=f)
    
    return(m)
}

read.smap = function(f=NULL) {
# This function will read in the cmap file as a smap data frame

    m = read.bionano.map(f=f)
    
    return(m)
}

read.map = function(f=NULL) {
# This function will read in the map file, which is not officially supported.
# map file format is depreciated but sometime still have useful information about alignment

    # first we read in header lines started with #
    header.obj = read.bionano.header(f)
    if (!is.null(header.obj$error)) { 
        result$error = -100
        result$errorMsg = paste("Unable to read in file '", f , "'", sep="")
        return(result)
    } 
    
    # we will read in map as data frame, since this is depreciated, I use hard coded line to get header.
	m = as.data.frame(read.table(f, header=T, skip = header.obj$headerLineCount, fill=T))
    
    return(m)
}

write.cmap = function(file=NULL, cmap.matrix = NULL) {
    # This is hard coded to write simple version of cmap without extra columns
    header1="# Generated by RScript write.cmap
# CMAP File Version:    0.1
# Label Channels:       1
# Nickase Recognition Site 1:   unknown
# Number of Consensus Nanomaps:    "
    header2="#h\tCMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence
#f\tint\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint"
    header = paste(header1, length(unique(cmap.matrix[,1])), "\n", header2, sep="")
    write(header, file=file)
    write.table(cmap.matrix, file=file, append=TRUE, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
}

write.xmap = function(file=NULL, xmap.matrix = NULL, header = NULL) {
    if (is.null(header)) {
    header="# Generated by RScript write.xmap
# XMAP File Version:	0.1
# Reference Maps From:	/mnt/data_out/human_ips_test_ips1_run1/contigs/exp_refineFinal1/exp_refineFinal1_contig1_r.cmap
# Query Maps From:	/mnt/data_out/human_ips_test_ips1_run1/contigs/exp_refineFinal1/exp_refineFinal1_contig1_q.cmap
#h XmapEntryID	QryContigID	RefcontigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum
#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 
"
    }
    write(header, file=file)
    write.table(xmap.matrix, file=file, append=TRUE, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
}


# Extended methods, internal use only
# Extract sub matrix from map matrix, used to generate label alignment
get.map.align.labID.matrix = function(map.m = NULL, col.prefix = "NanoLabelID") {
    # by default, this will extract NanoLabelID1.2.3....N as a matrix, if col.prefix is set to "RefLabelIDleft" or "RefLabelIDright", it will get these columns instead
    # this is used to align labels 
    # get.map.align.labID.matrix(map.m = all.ngs2bionano.map.m, col.prefix = "NanoLabelID")
    
    # First we find the last column and name and figure out how many 3X columns there are, RefLabelIDRight4 means that max is 4
    last.col.name = colnames(map.m)[NCOL(map.m)]
    n.cols = as.numeric(gsub("(\\D)", "", last.col.name, perl=TRUE))
    
    #colnames(map.m)[(NCOL(map.m)-3*n.cols) + 3*(1:n.cols)]
    selected.col.idx = grep(paste(col.prefix, "\\d+", sep=""), colnames(map.m))
    if (length(selected.col.idx) != n.cols) { stop("Unable to match regular expression of map.m colnames with expeteced pattern!") }
    
    return(map.m[,selected.col.idx])
}


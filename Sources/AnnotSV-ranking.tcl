############################################################################################################
# AnnotSV 2.1                                                                                              #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2019 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
#                                                                                                          #
# This is part of AnnotSV source code.                                                                     #
#                                                                                                          #
# This program is free software; you can redistribute it and/or                                            #
# modify it under the terms of the GNU General Public License                                              #
# as published by the Free Software Foundation; either version 3                                           #
# of the License, or (at your option) any later version.                                                   #
#                                                                                                          #
# This program is distributed in the hope that it will be useful,                                          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


## Define the position of some columns in an AnnotSV output line (columns requested for the ranking) 
proc SVprepareRanking {L_header} {

    global g_AnnotSV
    global g_i
    global L_Candidates

    # List of the candidate genes (given by the user)
    set L_Candidates {}
    if {$g_AnnotSV(candidateGenesFile) ne ""} {
	set L_Candidates [split [ContentFromFile $g_AnnotSV(candidateGenesFile)] " |\n|\t"]
    }

    # Check if we have all the needed information for ranking
    # Before this checking, ranking can not be done
    set g_AnnotSV(ranking) 0

    set Ls [split $L_header "\t"]

    # The user is using a VCF SV input file: the svtBEDcol is necessarily the 6th in the corresponding created annotated file
    if {[regexp ".vcf$|.vcf.gz$" $g_AnnotSV(SVinputFile)]} {
	set g_AnnotSV(svtBEDcol) 5
    } 

    # Check if we have all the needed information for ranking
    set g_i(gene)     [lsearch -regexp $Ls "Gene name"];    if {$g_i(gene) == -1} {puts "Gene name column not found => no SV ranking"; unset g_i; return}          
    set g_i(full)     [lsearch -regexp $Ls "AnnotSV type"]; if {$g_i(full) == -1} {puts "AnnotSV type column not found => no SV ranking"; unset g_i; return}        

    set g_i(GAINtot)  [lsearch -regexp $Ls "DGV_GAIN_n_samples_tested"]; if {$g_i(GAINtot) == -1} {puts "DGV_GAIN_n_samples_tested column not found => no SV ranking"; unset g_i; return}  
    set g_i(GAINfreq) [lsearch -regexp $Ls "DGV_GAIN_Frequency"];        if {$g_i(GAINfreq) == -1} {puts "DGV_GAIN_Frequency column not found => no SV ranking"; unset g_i; return}  
    set g_i(LOSStot)  [lsearch -regexp $Ls "DGV_LOSS_n_samples_tested"]; if {$g_i(LOSStot) == -1} {puts "DGV_LOSS_n_samples_tested column not found => no SV ranking"; unset g_i; return}  
    set g_i(LOSSfreq) [lsearch -regexp $Ls "DGV_LOSS_Frequency"];        if {$g_i(LOSSfreq) == -1} {puts "DGV_LOSS_Frequency column not found => no SV ranking"; unset g_i; return}  

    set g_i(dbVar_event)   [lsearch -regexp $Ls "dbVar_event"];   if {$g_i(dbVar_event) == -1} {puts "dbVar_event column not found => no SV ranking"; unset g_i; return}  
    set g_i(dbVar_status)  [lsearch -regexp $Ls "dbVar_status"];  if {$g_i(dbVar_status) == -1} {puts "dbVar_status column not found => no SV ranking"; unset g_i; return}  
    
    set g_i(morbidGenes)           [lsearch -regexp $Ls "morbidGenes"];           if {$g_i(morbidGenes) == -1} {puts "morbidGenes column not found => no SV ranking"; unset g_i; return}  
    set g_i(morbidGenesCandidates) [lsearch -regexp $Ls "morbidGenesCandidates"]; if {$g_i(morbidGenesCandidates) == -1} {puts "morbidGenesCandidates column not found => no SV ranking"; unset g_i; return}  

    set g_i(GHgene_elite)          [lsearch -regexp $Ls "GHgene_elite"];         
    set g_i(GHgene_not_elite)      [lsearch -regexp $Ls "GHgene_not_elite"];     

    set g_i(pLI)       [lsearch -regexp $Ls "pLI"];       if {$g_i(pLI) == -1} {puts "pLI column not found => no SV ranking"; unset g_i; return}  
    set g_i(HI_CGscore)   [lsearch -regexp $Ls "HI_CGscore"];   if {$g_i(HI_CGscore) == -1} {puts "HI_CGscore column not found => no SV ranking"; unset g_i; return}  
    set g_i(TriS_CGscore) [lsearch -regexp $Ls "TriS_CGscore"]; if {$g_i(TriS_CGscore) == -1} {puts "TriS_CGscore column not found => no SV ranking"; unset g_i; return}  
    
    # If we have all the needed information, ranking will be done
    set g_AnnotSV(ranking) 1
}


# If enhancer annotations is available, check the gene-enhancer relations (look at the gene, if it is in a precise list)
# Return:
#########
# - MorbidGenes            : SV overlaps an enhancer associated to a morbid gene                             (<=> ranking = 4)
# - DEL                    : SV overlaps an enhancer of a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2 (<=> ranking = 4)
# - DUP                    : SV overlaps the enhancer of a gene with a TriS_CGscore of 3 or 2                (<=> ranking = 4)
# - MorbidGenesCandidates  : SV overlaps an enhancer associated to a morbid gene candidate                   (<=> ranking = 3)
# - Candidates             : SV overlaps an enhancer associated to a gene candidate (given by the user)      (<=> ranking = 3)
# - NotUsed                : No GeneHancer annotation available
proc EnhancerInformation {Ls SVtype SVtoAnn} {
    global g_AnnotSV
    global g_i

    global rankingPreparation
    global L_MorbidGenes
    global L_DEL
    global L_DUP
    global L_MorbidGenesCandidates
    global L_Candidates

    global EliteGene
    global NotEliteGene


    # Do we have enhancer information?
    ##################################
    if {$g_i(GHgene_elite) == -1 || $g_i(GHgene_not_elite) == -1} {
	return "NotUsed"
    }

    # Creation of different genes list to prepare the ranking:
    ##########################################################
    if {![info exists rankingPreparation]} {
	set rankingPreparation 1

	# List of the morbid genes:
	regsub "Sources/?" $g_AnnotSV(sourcesDir) "Annotations/Genes-based/OMIM" omimDir
	set MorbidGenesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenes.tsv.gz"]
	set L_MorbidGenes {}
	if {[file exists $MorbidGenesFileFormattedGzip]} {
	    foreach L [LinesFromGZFile $MorbidGenesFileFormattedGzip] {
		lappend L_MorbidGenes [lindex $L 0]
	    }
	}
	# L_DEL: List of genes with a pLI > 0.9 or with HI_CGscore of 3 or 2
	# L_DUP: List of genes with a TriS_CGscore of 3 or 2
	regsub "Sources/?" $g_AnnotSV(sourcesDir) "Annotations/Genes-based/ExAC" ExACdir
	set ExACfile [glob -nocomplain "$ExACdir/*_GeneIntolerance.pLI-Zscore.annotations.tsv.gz"]
	set L_DEL {}
	if {[file exists $ExACfile]} {
	    foreach L [LinesFromGZFile $ExACfile] {
		set Ls [split $L "\t"]
		set pLI [lindex $Ls 3] 
		if {$pLI > 0.9} {
		    lappend L_DEL [lindex $Ls 0]
		}
	    }
	}
	regsub "Sources/?" $g_AnnotSV(sourcesDir) "Annotations/Genes-based/ClinGen" ClinGenDir
	set CGfile [glob -nocomplain "$ClinGenDir/*_ClinGenAnnotations.tsv.gz"]
	set L_DUP {}
	if {[file exists $CGfile]} {
	    foreach L [LinesFromGZFile $CGfile] {
		set Ls [split $L "\t"]
		set HI_CGscore [lindex $Ls 1] 
		if {$HI_CGscore eq "2" || $HI_CGscore eq "3"} {
		    lappend L_DEL [lindex $Ls 0]
		}
		set TriS_CGscore [lindex $Ls 2] 
		if {$TriS_CGscore eq "2" || $TriS_CGscore eq "3"} {
		    lappend L_DUP [lindex $Ls 0]
		}
 
	    }
	}
	# List of the morbid genes candidates:
	set MorbidGenesCandidatesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenescandidates.tsv.gz"]
	set L_MorbidGenesCandidates {}
	if {[file exists $MorbidGenesCandidatesFileFormattedGzip]} {
	    foreach L [LinesFromGZFile $MorbidGenesCandidatesFileFormattedGzip] {
		lappend L_MorbidGenesCandidates [lindex $L 0]
	    }
	}
    }


    # Listing of the elite_genes and not_elite_genes overlapped by the SV
    #####################################################################
    set L_enhancersAssociatedGenes {}
    set EG "[lindex $Ls $g_i(GHgene_elite)]"
    if {[regexp "\\.\\.\\.$" $EG] && [info exists EliteGene($SVtoAnn)]} {       ;# AnnotSV restrict the number of overlapping reported features to 20. Keep back all the genes values
	set EG $EliteGene($SVtoAnn)                                             ;# ok, there is no ";"
    } else {
	set EG [split $EG ";"]
    }
    set NEG "[lindex $Ls $g_i(GHgene_not_elite)]"
    if {[regexp "\\.\\.\\.$" $NEG] && [info exists NotEliteGene($SVtoAnn)]} {   ;# AnnotSV restrict the number of overlapping reported features to 20. Keep back all the genes values
	set NEG $NotEliteGene($SVtoAnn)                                         ;# ok, there is no ";"
    } else {
	set NEG [split $NEG ";"]
    } 
    set L_enhancersAssociatedGenes [lsort -unique [list {*}$EG {*}$NEG]]


    # Check if the SV overlaps an enhancer associated to a morbid gene 
    ##################################################################
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_MorbidGenes $g] ne -1} {
	    return "MorbidGenes"
	}
    }	


    # Check if the SV overlaps an enhancer associated to a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
    #########################################################################################################
    # Check for a del:
    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	foreach g $L_enhancersAssociatedGenes {
	    if {[lsearch -exact $L_DEL $g] ne -1} {
		return "DEL"
	    }
	}	    
    }

    # Check if the SV overlaps an enhancer associated to a gene with a TriS_CGscore of 3 or 2
    #########################################################################################
    # Check for a dup:
    if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	foreach g $L_enhancersAssociatedGenes {
	    if {[lsearch -exact $L_DUP $g] ne -1} {
		return "DUP"
	    }
	}
    }	    

    # Check if the SV overlaps an enhancer associated to a morbid gene candidate
    ############################################################################
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_MorbidGenesCandidates $g] ne -1} {
	    return "MorbidGenesCandidates"
	}
    }	

    # Check if the SV overlaps an enhancer associated to a candidate gene (given by the user)
    #########################################################################################
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_Candidates $g] ne -1} {
	    return "Candidates"
	}
    }	
}


## Rank the pathogenicity of the different SV as follows:
#########################################################
## category 1 = benign
##           > 70% SV overlapped with a benign SV + doesn't overlap with i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
## category 2 = likely benign
##           < 70% SV overlapped with a benign SV + doesn't overlap with i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene (or their enhancers)
## category 3 = VOUS (variant of unknown significance)
##           SV that overlaps i) a morbid gene candidate  (or its enhancer) and ii) a candidate gene  (or its enhancer) (with at least 1bp)
## category 4 = likely pathogenic 
##           SV that overlaps a morbid gene (or its enhancer) (with at least 1bp)
##           or for a del: SV that overlap a gene (or its enhancer) with a pLI > 0.9 or with HI_CGscore of 3 or 2
##           or for a dup: SV that overlap a gene (or its enhancer) TriS_CGscore of 3 or 2
## category 5 = pathogenic
##           SV that overlaps a pathogenic SV (with at least 1bp)

## ClinGen HI_CGscore and TriS_CGscore explanations:
####################################################
##   Rating	Possible Clinical Interpretation
##   ------     --------------------------------
##   3	        Sufficient evidence for dosage pathogenicity
##   2	        Some evidence for dosage pathogenicity
##   1	        Little evidence for dosage pathogenicity
##   0	        No evidence for dosage pathogenicity
##   40         Evidence suggests the gene is not dosage sensitive
##
## HI = Haploinsufficiency  TriS = Triplosensitivity
proc SVranking {L_annotations} {

    global g_AnnotSV
    global g_i
    global L_Candidates

    # Check if we have enougth information to do the ranking:
    #########################################################
    set ranking ""
    if {$g_AnnotSV(svtBEDcol) eq -1} {return $ranking}
    if {!$g_AnnotSV(ranking)} {return $ranking}


    # Ranking!!
    ###########

    set Ls [split $L_annotations "\t"]
    set SVtype [lindex $Ls $g_AnnotSV(svtBEDcol)]   
    set SVtoAnn [join [lrange $Ls 1 3] ","]
    set enhancer [EnhancerInformation $Ls $SVtype $SVtoAnn]

    ## category 5 = pathogenic
    ##              SV that overlap a pathogenic SV (with at least 1bp)
    ###################################################################
    set dbVar_status [lindex $Ls $g_i(dbVar_status)]
    if {$dbVar_status ne ""} {
	set dbVar_event [lindex $Ls $g_i(dbVar_event)]
	if {[regexp -nocase "Del|Loss|<CN0>" $SVtype] && [regexp -nocase "Del|Loss|<CN0>" $dbVar_event]} {
	    set ranking "5"	
	    return $ranking
	}
	if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype] && [regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $dbVar_event]} {
	    set ranking "5"	
	    return $ranking
	}
    }

    ## category 4 = likely pathogenic 
    ##           SV that overlaps a morbid gene (or its enhancer) (with at least 1bp)
    ##           or for a del: SV that overlaps a gene (or its enhancer) with a pLI > 0.9 or with HI_CGscore of 3 or 2
    ##           or for a dup: SV that overlaps a gene (or its enhancer) with a TriS_CGscore of 3 or 2
    ###################################################################

    # Check if a SV overlaps a morbid gene 
    set morbidGenes [lindex $Ls $g_i(morbidGenes)]
    if {[regexp "yes" $morbidGenes]} {
	set ranking "4"	
	return $ranking
    }
    # Check if a SV overlaps an enhancer associated to a morbid gene 
    if {$enhancer == "MorbidGenes"} {
	set ranking "4"	
	return $ranking
    }
    
    # Check for a del:
    set pLI [lindex $Ls $g_i(pLI)]
    set HI_CGscore [lindex $Ls $g_i(HI_CGscore)]
    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	# Check SV that overlap a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
	if {$pLI > 0.9 || $HI_CGscore eq 3 || $HI_CGscore eq 2} {    ; # {"" > 0.9} is false; code ok
	    set ranking "4"	
	    return $ranking
	}
	# Check SV that overlap the enhancer of a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
	if {$enhancer == "DEL"} {
	    set ranking "4"	
	    return $ranking
	}
    }
    
    # Check for a dup:
    set TriS_CGscore [lindex $Ls $g_i(TriS_CGscore)]
    if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	# Check SV that overlap a gene TriS_CGscore of 3 or 2
	if {$TriS_CGscore eq 3 || $TriS_CGscore eq 2} {
	    set ranking "4"	
	    return $ranking
	}
	# Check SV that overlap the enhancer of a gene TriS_CGscore of 3 or 2
	if {$enhancer == "DUP"} {
	    set ranking "4"	
	    return $ranking
	}
    }
    
    
    ## category 3 = VOUS (variant of unknown significance)
    ##           SV that overlap an enhancer or a CDS from i) a morbid gene candidate and ii) a candidate gene (with at least 1bp)
    ###################################################################

    # Check if a SV overlaps a morbid gene candidate
    set morbidGenesCandidates [lindex $Ls $g_i(morbidGenesCandidates)]
    if {[regexp "yes" $morbidGenesCandidates]} {
	set ranking "3"	
	return $ranking
    }
    # Check if a SV overlaps an enhancer associated to a morbid gene candidate
    if {$enhancer == "MorbidGenesCandidates"} {
	set ranking "3"	
	return $ranking
    }

    if {$g_AnnotSV(candidateGenesFile) ne ""} {
	# Check if a SV overlaps a CDS from a candidate gene (with at least 1bp)
	foreach gene [split [lindex $Ls $g_i(gene)] "/"] {
	    if {$gene ne "" && [lsearch -exact $L_Candidates $gene] ne -1} {
		set ranking "3"	
		return $ranking
	    }
	}
	# Check if a SV overlaps the enhancer from a candidate gene (with at least 1bp)
	if {$enhancer == "Candidates"} {
	    set ranking "3"	
	    return $ranking
	}
    }


    ## category 1 = benign
    ##           > 80% SV overlapped with a benign SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
    ###################################################################
    set GAINtot [lindex $Ls $g_i(GAINtot)] 
    regsub ","  [lindex $Ls $g_i(GAINfreq)] "." GAINfreq; # needed with -metrics=fr 
    set LOSStot [lindex $Ls $g_i(LOSStot)]
    regsub ","  [lindex $Ls $g_i(LOSSfreq)] "." LOSSfreq; # needed with -metrics=fr 
    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	if {$LOSStot > $g_AnnotSV(minTotalNumber) && $LOSSfreq > 0.01} {
	    set ranking "1"	
	    return $ranking
	}
    } elseif {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	if {$GAINtot > $g_AnnotSV(minTotalNumber) && $GAINfreq > 0.01} {
	    set ranking "1"	
	    return $ranking
	}
    }
    
    ## category 2 = likely benign
    ##           < 80% SV overlapped with a benign SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
    ###################################################################    
    set ranking "2"	
    return $ranking
}

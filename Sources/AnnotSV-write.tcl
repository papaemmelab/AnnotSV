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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################



proc OrganizeAnnotation {} {

    global g_AnnotSV
    global g_Lgenes
    global VCFheader
    global g_numberOfAnnotationCol
    global headerFileToRemove 
    global L_Candidates
    global g_SVLEN


    # OUTPUT
    ###############
    set tmpFullAndSplitBedFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp"
    set outputFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)" 


    # Check the -svtBEDcol option
    # Incr -svtBEDcol +1 (-1 to be ok for an informatic list + 2 for the 2 added columns)
    #####################################################################################
    checksvtBEDcol $g_AnnotSV(bedFile)


    ################### Display of the annotations to realize ########################
    ##################################################################################
    puts "\n...listing of the annotations to realized ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
    puts "\t...refGene annotation" 
    puts "\t(with $g_AnnotSV(refGene))"
    ####### "Genes-based annotations"
    puts "\t...Genes-based annotations"
    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} { ; # g_AnnotSV(extann) defined in the main
	ExternalAnnotations
	foreach F [ExternalAnnotations L_Files] {
	    # we remove the first "genes" column (not reported in output, used as an ID)
	    regsub "\[^\t\]+\t" [ExternalAnnotations $F Header] "" extannHeader
	    append allGenesBasedHeader "\t$extannHeader"
	}
    }
    ####### "SVincludedInFt"
    puts "\t...Annotations with features overlapping the SV"
    if {$g_AnnotSV(dgvAnn)} {puts "\t\t...DGV Gold Standard frequency annotation"}
    if {$g_AnnotSV(gnomADann)} {puts "\t\t...gnomAD frequency annotation"}
    if {$g_AnnotSV(DDDfreqAnn)} {puts "\t\t...DDD frequency annotation"}
    if {$g_AnnotSV(1000gAnn)} {puts "\t\t...1000g frequency annotation"}
    if {$g_AnnotSV(IMHann)} {puts "\t\t...Ira M. Hall's lab frequency annotation"}
    regsub "Sources/?" $g_AnnotSV(sourcesDir) "Annotations/Users/$g_AnnotSV(genomeBuild)" usersDir
    foreach formattedUserBEDfile [glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] {
	puts "\t\t...[file tail $formattedUserBEDfile]"
	regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	set L1header [split [FirstLineFromFile $userHeaderFile] "\t"]
	set L_headerColName [lrange $L1header 3 end]
	append SVincludedInFtHeader "\t[join $L_headerColName "\t"]"
	# Number of columns of annotation in the user bedfile (without the 3 columns "chrom start end")
	set nColHeader [llength $L_headerColName]
	puts "\t\t($nColHeader annotations columns: [join $L_headerColName ", "])"
    }
    ####### "FtIncludedInSV"
    puts "\t...Annotations with features overlapped with the SV"
    puts "\t\t...Promoters annotation"
    if {$g_AnnotSV(NRSVann)} {puts "\t\t...dbVar_pathogenic_NR_SV annotation"}
    if {$g_AnnotSV(GHann)} {puts "\t\t...GH annotation"}
    if {$g_AnnotSV(tadAnn)} {puts "\t\t...TAD annotation"}
    foreach formattedUserBEDfile [glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] {
	puts "\t\t...[file tail $formattedUserBEDfile]"
	regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	set L1header   [split [FirstLineFromFile $userHeaderFile] "\t"]
	set L_headerColName [lrange $L1header 3 end]
	append FtIncludedInSVHeader "\t[join $L_headerColName "\t"]"
	# Number of columns of annotation in the user bedfile (without the 3 columns "chrom start end")
	set nColHeader [llength $L_headerColName]
	puts "\t\t($nColHeader annotations columns: [join $L_headerColName ", "])"
    }
    ####### "Breakpoints annotations"
    puts "\t...Breakpoints annotations"
    if {$g_AnnotSV(gcContentAnn)} {puts "\t\t...GC content annotation"}
    if {$g_AnnotSV(repeatAnn)} {puts "\t\t...Repeat annotation"}




    ################### Display: annotation in progress ####################
    ########################################################################
    puts "\n\n...annotation in progress ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 



    ################### Writing of the header (first line of the output) ####################
    #########################################################################################
    set headerOutput "AnnotSV ID\tSV chrom\tSV start\tSV end\tSV length\t"
    if {[info exist VCFheader]} {
	# SVinputFile = VCF
	append headerOutput "$VCFheader\t"
	set theLength [expr {[llength [split $headerOutput "\t"]]-3}] ; # We remove "1" for the count of the last "\t" + we remove "2" for the 2 columns not in the input file (AnnotSV ID and SV length)
	set g_AnnotSV(SVinputInfo) 1 ; # The created bedfile contains only the info to report
    } else {
	# SVinputFile = BED
	if {$g_AnnotSV(SVinputInfo)} {
	    # The user wants to keep all the columns from the SV BED input file
	    regsub -nocase ".bed$" $g_AnnotSV(SVinputFile) ".header.tsv" headerSVinputFile
	    set theLength [llength [split [FirstLineFromFile $g_AnnotSV(bedFile)] "\t"]] ;# theLength = number of col in the SV input BED file
	    set testH 1
	    if {[file exists $headerSVinputFile]} {
		# The user has given a header for all the columns from the SV BED input file
		set headerFromTheUser [split [FirstLineFromFile $headerSVinputFile] "\t"]	
		if {[llength $headerFromTheUser] ne $theLength} {
		    puts "Numbers of columns from $g_AnnotSV(SVinputFile) and $headerSVinputFile are different ($theLength != [llength $headerFromTheUser])" 
		    puts "=> Can not report: $headerFromTheUser"
		} else {
		    if {$g_AnnotSV(svtBEDcol) ne -1} {
			set i [expr {$g_AnnotSV(svtBEDcol)-2}]
			set headerFromTheUser [lreplace $headerFromTheUser $i $i "SV type"]
		    }
		    append headerOutput "[join [lrange $headerFromTheUser 3 end] "\t"]\t"		
		    set testH 0
		}	
	    } 
	    if {$testH} {
		set i 5 ; #AnnotSV ID	SV chrom    SV start	SV end	SV length
		set j [expr {$theLength+2}]
		while {$i < $j} {
		    if {$i eq $g_AnnotSV(svtBEDcol)} {
			append headerOutput "SV type\t"
		    } else {
			append headerOutput "\t"
		    }
		    incr i
		}
	    }
	} elseif {$g_AnnotSV(svtBEDcol) ne -1} {; # At least the "SV type" column should be reported for the ranking
	    # The user doesn't want to keep all the columns from the SV BED input file. We keep only the SV type (for the ranking)
	    append headerOutput "SV type\t"
	}
    }

    append headerOutput "AnnotSV type\tGene name\tNM\tCDS length\ttx length\tlocation\tintersectStart\tintersectEnd"

    ####### "SVincludedInFt"
    if {$g_AnnotSV(dgvAnn)} {
	append headerOutput "\tDGV_GAIN_IDs\tDGV_GAIN_n_samples_with_SV\tDGV_GAIN_n_samples_tested\tDGV_GAIN_Frequency"
	append headerOutput "\tDGV_LOSS_IDs\tDGV_LOSS_n_samples_with_SV\tDGV_LOSS_n_samples_tested\tDGV_LOSS_Frequency"
    }
    if {$g_AnnotSV(gnomADann)} {
	append headerOutput "\tGD_ID\tGD_AN\tGD_N_HET\tGD_N_HOMALT\tGD_AF\tGD_POPMAX_AF\tGD_ID_others"
    }
    if {$g_AnnotSV(DDDfreqAnn)} {
	append headerOutput "\tDDD_SV\tDDD_DUP_n_samples_with_SV\tDDD_DUP_Frequency\tDDD_DEL_n_samples_with_SV\tDDD_DEL_Frequency"
    }
    if {$g_AnnotSV(1000gAnn)} {
	append headerOutput "\t1000g_event\t1000g_AF\t1000g_max_AF"
    }
    if {$g_AnnotSV(IMHann)} {
	append headerOutput "\tIMH_ID\tIMH_AF\tIMH_ID_others"
    }
    if {[glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] ne ""} {
	append headerOutput "$SVincludedInFtHeader"
    }

    ####### "FtIncludedInSV"
    append headerOutput "\tpromoters"
    if {$g_AnnotSV(NRSVann)} {
	append headerOutput "\tdbVar_event\tdbVar_variant\tdbVar_status"
    }
    if {$g_AnnotSV(GHann)} {
	append headerOutput "\tGHid_elite\tGHid_not_elite\tGHtype\tGHgene_elite\tGHgene_not_elite\tGHtissue"
    }
    if {$g_AnnotSV(tadAnn)} {
	append headerOutput "\tTADcoordinates\tENCODEexperiments"
    }
    if {$g_AnnotSV(vcfFiles) ne ""} {
	foreach sample $g_AnnotSV(vcfSamples) {
	    append headerOutput "\t#hom($sample)\t#htz($sample)"
	}
    }
    if {$g_AnnotSV(filteredVCFfiles) ne ""} {
	foreach sample $g_AnnotSV(filteredVCFsamples) {
	    append headerOutput "\tcompound-htz($sample)"
	}
    }
    if {[glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] ne ""} {
	append headerOutput "$FtIncludedInSVHeader"
    }

    ####### "Breakpoints annotations"
    if {$g_AnnotSV(gcContentAnn)} {
	append headerOutput "\tGCcontent_left\tGCcontent_right"
    }
    if {$g_AnnotSV(repeatAnn)} {
	append headerOutput "\tRepeats_coord_left\tRepeats_type_left\tRepeats_coord_right\tRepeats_type_right"
    }

    ####### "Genes-based annotations"
    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} { ; # g_AnnotSV(extann) defined in the main
	append headerOutput "$allGenesBasedHeader"
    }


    
    append headerOutput "\tAnnotSV ranking"
    ReplaceTextInFile $headerOutput $outputFile





    # Preparation for the ranking (from benign to pathogenic)
    #########################################################
    SVprepareRanking $headerOutput    ; # svtBEDcol (for VCF input file) is defined there

    if {$g_AnnotSV(svtBEDcol) eq -1} { ; # SV type is required for the ranking
	puts "\nWARNING: AnnotSV requires the SV type (duplication, deletion...) to classify the SV"
	puts "Not provided (svtBEDcol = -1)"
	puts "=> No SV ranking"
    }


    # Number of columns from each Genes-based file
    ##############################################
    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} {		    
	foreach F [ExternalAnnotations L_Files] {
	    set nbColumns($F) [expr {[llength [split [ExternalAnnotations $F Header] "\t"]]-1}]
	}
    }
   

    # Prepare the intersection between SV and BED file (DGV, dbVar...):
    # => creation of a bedfile (chrom, start and end coordinates) with the SV to annotate + the intersections with genes
    ####################################################################################################################
    set g_AnnotSV(fullAndSplitBedFile) "$g_AnnotSV(outputDir)/[file tail $g_AnnotSV(bedFile)].users.bed"
    file delete -force $g_AnnotSV(fullAndSplitBedFile)
    set L_UsersText {}
    foreach L [LinesFromFile $tmpFullAndSplitBedFile] {
	set Ls [split $L "\t"]
	set AnnotSVtype [lindex $Ls end]
	if {$AnnotSVtype eq "split"} {
	    set SVchrom [lindex $Ls 0]
	    set SVleft  [lindex $Ls 1]
	    set SVright [lindex $Ls 2]
	    set exonStarts [lindex $Ls end-4]
	    set exonEnds   [lindex $Ls end-3]
	    set tx_left [lindex [split $exonStarts ","] 0]
	    set tx_right [lindex [split $exonEnds ","] end-1]
	    if {$SVleft<$tx_left} {set intersectStart "$tx_left"} else {set intersectStart "$SVleft"}
	    if {$SVright<$tx_right} {set intersectEnd "$SVright"} else {set intersectEnd "$tx_right"}
	    lappend L_UsersText "$SVchrom\t$intersectStart\t$intersectEnd"
	} else {
	    lappend L_UsersText "[join [lrange $Ls 0 2] "\t"]"
	}
    }
    set L_UsersText [lsort -unique $L_UsersText]
    WriteTextInFile [join $L_UsersText "\n"] $g_AnnotSV(fullAndSplitBedFile)
    checkBed $g_AnnotSV(fullAndSplitBedFile)
    regsub -nocase ".bed$" $g_AnnotSV(fullAndSplitBedFile) ".formatted.bed" newFullAndSplitBedFile

    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    set g_AnnotSV(fullAndSplitBedFile) "$g_AnnotSV(outputDir)/[file tail $g_AnnotSV(bedFile)].users.sorted.bed"
    if {[catch {eval exec sort -k1,1 -k2,2n $newFullAndSplitBedFile > $g_AnnotSV(fullAndSplitBedFile)} Message]} {
	puts "-- OrganizeAnnotation --"
	puts "sort -k1,1 -k2,2n $newFullAndSplitBedFile > $g_AnnotSV(fullAndSplitBedFile)"
	puts "$Message"
	puts "Exit with error"
	exit 2
    }

    file delete -force $newFullAndSplitBedFile
    unset L_UsersText


    # Parse
    ###############
    set L_TextToWrite {}
    foreach L [LinesFromFile $tmpFullAndSplitBedFile] {

	set Ls [split $L "\t"]

	# Full + split
	set SVchrom   [lindex $Ls 0]
	set SVleft    [lindex $Ls 1]
	set SVright   [lindex $Ls 2]
	set AnnotSVtype   [lindex $Ls end]                 ;# full or split
	if {$g_AnnotSV(svtBEDcol) ne -1} { 
	    set i [expr {"$g_AnnotSV(svtBEDcol)"-2}]
	    set SVtype [lindex $Ls $i]                     ;# DEL, DUP, <CN0>...
	} else {
	    set SVtype ""
	}
	
	# For the use of the -typeOfAnnotation option: keep only the corresponding lines (full or split or both)
	if {$g_AnnotSV(typeOfAnnotation) eq "full" && $AnnotSVtype eq "split"} {continue}
	if {$g_AnnotSV(typeOfAnnotation) eq "split" && $AnnotSVtype eq "full"} {continue}


	if {$AnnotSVtype eq "split"} {
	    # split
	    set txStart    [lindex $Ls end-11]
	    set txEnd      [lindex $Ls end-10]
	    set strand     [lindex $Ls end-9]
	    set geneName   [lindex $Ls end-8]
	    set NM         [lindex $Ls end-7]
	    set CDSstart   [lindex $Ls end-6]
	    set CDSend     [lindex $Ls end-5]
	    set exonStarts [lindex $Ls end-4]
	    set exonEnds   [lindex $Ls end-3]
	    set CDSl       [lindex $Ls end-2]
	    set txL        [lindex $Ls end-1]
	    set locationStart ""
	    set locationEnd ""
	    set intersectStart ""
	    set intersectEnd ""
	} else {
	    set SV "[join [lrange $Ls 0 2] "\t"]"
	    # full
	    if {[info exists g_Lgenes($SV)]} {
		set geneName "$g_Lgenes($SV)"  ; # SV/oldSV defined with: "chrom, start, end, SVtype":	    set NM         ""
	    } else {
		set geneName ""
	    }
	    set NM         ""
	    set CDSl       ""
	    set txL        ""
	    set location   ""
	    set intersect  "\t"
	    set compound   ""
	}	    

	if {$g_AnnotSV(candidateGenesFiltering) eq "yes"} {
	    set test 0
	    foreach g [split $geneName "/"] {
		if {[lsearch -exact $L_Candidates $g] eq -1} {
		    set test 1
		}
	    }	
	    if {$test} {continue}
	}

	# Definition of "locationStart" and "locationEnd" variables
	if {$AnnotSVtype eq "split"} {
	    set nbEx [expr {[llength [split $exonStarts ","]]-1}] ; #Example: 1652370,1657120,1660664,1661968 --> 1652370 1657120 1660664 1661968 {}
	    
	    set tx_left [lindex [split $exonStarts ","] 0]
	    set tx_right [lindex [split $exonEnds ","] end-1]
	    
	    if {$SVleft<=$tx_left} {         ; # SV begins before tx start (or tx end for strand "-")
		if {$strand eq "+"} {
		    set locationStart "txStart"
		} else {
		    set locationEnd "txEnd"
		}
	    }
	    set i 0
	    foreach A [split $exonStarts ","] B [split $exonEnds ","] {
		if {$A eq "" || $B eq ""} {continue}
		incr i
		
		# SV left
		if {$SVleft<$B} {
		    if {$SVleft>$A} {
			if {$strand eq "+"} {
			    set locationStart "exon$i"
			} else {
			    set locationEnd "exon[expr {$nbEx-$i+1}]" ; # gene on the strand "-"
			}
		    } else {
			if {$strand eq "+"} {
			    if {$locationStart eq ""} {set locationStart "intron[expr {$i-1}]"}
			} else {
			    if {$locationEnd eq ""} {set locationEnd "intron[expr {$nbEx-$i+1}]"}
			}
		    }
		}
		# SV right	
		if {$SVright<$B} {
		    if {$SVright>$A} {
			if {$strand eq "+"} {
			    set locationEnd "exon$i"; break
			} else {
			    set locationStart "exon[expr {$nbEx-$i+1}]"; break
			}
		    } else {
			if {$strand eq "+"} {
			    set locationEnd "intron[expr {$i-1}]"; break
			} else {
			    set locationStart "intron[expr {$nbEx-$i+1}]"; break
			}
		    }
		}
	    }
	    if {$locationEnd eq "" || $locationStart eq ""} {     ; # SV finishes after tx end
		if {$strand eq "+"} {
		    set locationEnd "txEnd"
		} else {
		    set locationStart "txStart"
		}
	    }
	    set location "$locationStart-$locationEnd"
	} 

	# Definition of "intersectStart" and "intersectEnd" variables
	if {$AnnotSVtype eq "split"} {
	    if {$SVleft<$tx_left} {set intersectStart "$tx_left"} else {set intersectStart "$SVleft"}
	    if {$SVright<$tx_right} {set intersectEnd "$SVright"} else {set intersectEnd "$tx_right"}
	    set intersect "$intersectStart\t$intersectEnd"
	} 

	# Promoters annotation
	set promoterText ""
	if {$g_AnnotSV(promAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set promoterText "[promoterAnnotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set promoterText "[promoterAnnotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# DGV annotation
	set dgvText ""
	if {$g_AnnotSV(dgvAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set dgvText "[DGVannotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set dgvText "[DGVannotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# gnomAD annotations
	set gnomADtext ""
	if {$g_AnnotSV(gnomADann)} {
	    if {$AnnotSVtype eq "split"} {
		set gnomADtext "[gnomADannotation $SVchrom $intersectStart $intersectEnd $SVtype]"
	    } else {
		set gnomADtext "[gnomADannotation $SVchrom $SVleft $SVright $SVtype]"
	    } 
	}

	# DDD frequency + gene annotations
	set dddText ""
	if {$g_AnnotSV(DDDfreqAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set dddText "[DDDfrequencyAnnotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set dddText "[DDDfrequencyAnnotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# 1000g annotations
	set 1000gText ""
	if {$g_AnnotSV(1000gAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set 1000gText "[1000gAnnotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set 1000gText "[1000gAnnotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# Ira M. Hall's lab annotations
	set IMHtext ""
	if {$g_AnnotSV(IMHann)} {
	    if {$AnnotSVtype eq "split"} {
		set IMHtext "[IMHannotation $SVchrom $intersectStart $intersectEnd $SVtype]"
	    } else {
		set IMHtext "[IMHannotation $SVchrom $SVleft $SVright $SVtype]"
	    } 
	}

	# dbVar pathogenic NR SV annotation
	set NRSVtext ""
	if {$g_AnnotSV(NRSVann)} {
	    if {$AnnotSVtype eq "split"} {
		set NRSVtext "[pathogenicNRSVannotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set NRSVtext "[pathogenicNRSVannotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# User FtIncludedInSV BED annotations. 
	set L_FtIncludedInSVtext {}
   	foreach formattedUserBEDfile [glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] {
	    if {$AnnotSVtype eq "split"} {
		lappend L_FtIncludedInSVtext "[userBEDannotation $formattedUserBEDfile $SVchrom $intersectStart $intersectEnd]"
	    } else {
		lappend L_FtIncludedInSVtext "[userBEDannotation $formattedUserBEDfile $SVchrom $SVleft $SVright]"
	    } 
	}
	set FtIncludedInSVtext [join $L_FtIncludedInSVtext "\t"]

	# User SVincludedInFt BED annotations. 
	set L_SVincludedInFtText {}
   	foreach formattedUserBEDfile [glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] {
	    if {$AnnotSVtype eq "split"} {
		lappend L_SVincludedInFtText "[userBEDannotation $formattedUserBEDfile $SVchrom $intersectStart $intersectEnd]"
	    } else {
		lappend L_SVincludedInFtText "[userBEDannotation $formattedUserBEDfile $SVchrom $SVleft $SVright]"
	    } 
	}
	set SVincludedInFTtext [join $L_SVincludedInFtText "\t"]

	# Genes-based annotations ($genesBasedText) in the output file
	set genesBasedText ""
	if {$AnnotSVtype eq "split"} {
	    ######### split lines ###################################
	    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) != ""} {		    
		foreach F [ExternalAnnotations L_Files] {
		    set AnnotFound "[ExternalAnnotations $F $geneName]"		
		    if {$AnnotFound eq ""} {
			lappend genesBasedText "[join [lrepeat $nbColumns($F) ""] "\t"]"
		    } else {		
			# Change metrics from "." to ","
			if {[set g_AnnotSV(metrics)] eq "fr"} {
			    foreach valueByColumn [split "$AnnotFound" "\t"] {
				if {[regexp "^(-)?\[0-9\]+\\.\[0-9\]+$" $valueByColumn]} { ;# Only for the values like "0.5", "-25.3" 
				    regsub -all {\.} $valueByColumn "," valueByColumn				    
				}
				lappend genesBasedText "$valueByColumn"
			    } 
			} else {
			    lappend genesBasedText "$AnnotFound"
			}
		    }
		}
		set genesBasedText [join $genesBasedText "\t"]
	    }
	} else {
	    ######### full lines ###################################
	    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) != ""} {	
		if {$geneName eq ""} {
		    foreach F [ExternalAnnotations L_Files] {
			lappend genesBasedText "[join [lrepeat $nbColumns($F) ""] "\t"]"
		    }
		    set genesBasedText [join $genesBasedText "\t"]
		} else {
		    set allGenesFromFullLine [split $geneName "/"]
		    foreach F [ExternalAnnotations L_Files] {
			
			# First, we search for the annotation of each gene that we merge with a "/"
			set L_AnnotFound {}
			foreach g $allGenesFromFullLine {
			    set AnnotFound "[ExternalAnnotations $F $g]"		
			    if {$AnnotFound eq ""} {
				set AnnotFound "[join [lrepeat $nbColumns($F) ""] "\t"]"
			    } 
			    lappend L_AnnotFound "$AnnotFound"		
			}
			set tutu [MergeAnnotation $L_AnnotFound $nbColumns($F)]
			# Second (for full lines), we doesn't keep the annotation of all genes (value is set to empty), 
			# except for scores and percentages where we keep the max value
			# (in order not to have value such "1.5/////-0.2////")
			set newGenesBasedText ""
			if {$tutu eq ""} {lappend genesBasedText ""} ; #if {$tutu eq ""} => doesn't enter in the foreach
			foreach valueByColumn [split $tutu "\t"] { 
			    if {$valueByColumn ne ""} {
				set isAScore 1
				foreach valueByGene [split $valueByColumn "/"] {
				    if {[regexp "ClinGenAnnotations.tsv" $F]} {
					set max ""
					# For ClinGen file (HI_CGscore + TriS_CGscore), the values are ordered as follow: 
					# 3 > 2 > 1 > 0 > 40 > 30 > Not yet evaluated
					# We only report 3, 2 or 1 in a full line
					if {[lsearch -exact {3 2 1} $valueByGene] ne -1} {
					    if {$valueByGene eq 3} {set max 3; break}
					    if {$valueByGene eq 2} {
						set max 2
					    } elseif {$max ne 2} {
						set max 1
					    } 
					}
				    } elseif {[regexp "morbid" $F]} {
					if {[regexp "yes" $valueByGene]} {set max "yes"}
				    } else {
					set max -1000
					if {[regexp "^(-)?\[0-9\]{1,3}\\.\[0-9\]+$" $valueByGene]} { ;# Only for the values like "0.25", "-25.3" 
					    if {$valueByGene > $max} {set max $valueByGene}
					} else {set isAScore 0; break}
				    } 
				}
				if {$isAScore} {
				    # Change metrics from "." to ","
				    if {[set g_AnnotSV(metrics)] eq "fr"} {
					regsub -all {\.} $max "," max
				    }
				    lappend newGenesBasedText $max
				} else {
				    lappend newGenesBasedText ""
				}
			    } else {
				lappend newGenesBasedText ""
			    }
			}
			lappend genesBasedText {*}$newGenesBasedText
		    }
		    set genesBasedText [join $genesBasedText "\t"]
		}
	    }
	}

	# GC content annotation
	set gcContentText ""
	if {$g_AnnotSV(gcContentAnn)} {
	    if {$AnnotSVtype ne "split"} {	
		set gcContentText "[GCcontentAnnotation $SVchrom $SVleft]"
		append gcContentText "\t[GCcontentAnnotation $SVchrom $SVright]"
	    } else {set gcContentText "\t"}
	}

	# Repeat annotation
	set repeatText ""
	if {$g_AnnotSV(repeatAnn)} {
	    if {$AnnotSVtype ne "split"} {	
		set repeatText "[RepeatAnnotation $SVchrom $SVleft]"
		append repeatText "\t[RepeatAnnotation $SVchrom $SVright]"
	    } else {set repeatText "\t\t\t"}
	}

	# GH annotation
	set GHtext ""
	if {$g_AnnotSV(GHann)} {
	    if {$AnnotSVtype eq "split"} {
		set GHtext "[GHannotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set GHtext "[GHannotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# TAD annotation
	set tadText ""
	if {$g_AnnotSV(tadAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set tadText "[TADannotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set tadText "[TADannotation $SVchrom $SVleft $SVright]"
	    } 
	}

	# Calculate the #hom and #htz variables
	set HomHtz ""
	if {$g_AnnotSV(vcfFiles) ne ""} {
	    if {$AnnotSVtype eq "split"} {
		set HomHtz "[VCFannotation $SVchrom $intersectStart $intersectEnd]"
	    } else {
		set HomHtz "[VCFannotation $SVchrom $SVleft $SVright]"
	    } 
	} 

	# Calculate the compound-htz variable
	if {$AnnotSVtype eq "split"} {
	    set compound ""
	    if {$g_AnnotSV(filteredVCFfiles) ne ""} {
		set compound "[filteredVCFannotation $SVchrom $txStart $txEnd]"
	    } 
	}

	# "bestAnn" annotation order: 
	# chrom txStart txEnd name2 name cdsStart cdsEnd exonStarts exonEnds
	#
	# headerOutput:
	#  "Gene name\tNM\tCDS length\ttx length\tlocation\tintersectStart\tintersectEnd\tOMIM ID\tOMIM phenotype\tOMIM inheritance\t#hom\t#htz"

	# Insertion of the SV length in the fourth column:
	set SVchrom [lindex $Ls 0]
	set SVstart [lindex $Ls 1]
	set SVend [lindex $Ls 2]
	# Creation of the AnnotSV ID (chrom_start_end_SVtype)
	set ID "${SVchrom}_${SVstart}_${SVend}_$SVtype"
	# Report of the SV length
	#set SVlength [expr {$SVend-$SVstart}] ; # No! Wrong for an insertion, a BND or a translocation
	if {[info exists g_SVLEN($ID)]} {
	    set SVlength $g_SVLEN($ID)
	} else {
	    if {[regexp "DEL" [normalizeSVtype $SVtype]]} { ;# DEL
		set SVlength [expr {$SVstart-$SVend}]
	    } elseif {[normalizeSVtype $SVtype] eq "INV" || $SVtype eq "DUP" || [regexp -nocase "<CN(\[0-9\]+)>" $SVtype]} { ;# DUP or INV
		set SVlength [expr {$SVend-$SVstart}]
	    } else {set SVlength ""}
	}

	################ Writing
	if {$g_AnnotSV(SVinputInfo)} {
	    set toadd [lrange $Ls 0 [expr {$theLength-1}]]
	    set toadd [linsert $toadd 3 $SVlength]
	    set TextToWrite "$ID\t[join $toadd "\t"]\t$AnnotSVtype\t$geneName\t$NM\t$CDSl\t$txL\t$location\t$intersect"
	} else {
	    if {$g_AnnotSV(svtBEDcol) ne -1} { ; # SV type is required for the ranking
		set TextToWrite "$ID\t[join [lrange $Ls 0 2] "\t"]\t$SVlength\t$SVtype\t$AnnotSVtype\t$geneName\t$NM\t$CDSl\t$txL\t$location\t$intersect"
	    } else {
		set TextToWrite "$ID\t[join [lrange $Ls 0 2] "\t"]\t$SVlength\t$AnnotSVtype\t$geneName\t$NM\t$CDSl\t$txL\t$location\t$intersect"
	    }
	}
	####### "SVincludedInFt"
	if {$dgvText ne ""} {
	    append TextToWrite "\t$dgvText"
	}
	if {$gnomADtext ne ""} {
	    append TextToWrite "\t$gnomADtext"
	}
	if {$dddText ne ""} {
	    append TextToWrite "\t$dddText"
	}
	if {$1000gText ne ""} {
	    append TextToWrite "\t$1000gText"
	}
	if {$IMHtext ne ""} {
	    append TextToWrite "\t$IMHtext"
	}
	if {[glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] ne ""} { ; # Don't put {$SVincludedInFTtext ne ""}: the user BED could have only 1 annotation column, and so $UserText can be equel to "" (without "\t")
	    append TextToWrite "\t$SVincludedInFTtext"
	}
	####### "FtIncludedInSV"
	append TextToWrite "\t$promoterText"
	if {$NRSVtext ne ""} {
	    append TextToWrite "\t$NRSVtext"
	}
	if {$GHtext ne ""} {
	    append TextToWrite "\t$GHtext"
	}
	if {$tadText ne ""} {
	    append TextToWrite "\t$tadText"
	}
	if {$HomHtz ne ""} {
	   append TextToWrite "\t$HomHtz"
	} 
	if {$compound ne ""} {
	   append TextToWrite "\t$compound"
	}
	if {[glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] ne ""} { ; # Don't put {$FtIncludedInSVtext ne ""}: the user BED could have only 1 annotation column, and so $UserText can be equel to "" (without "\t")
	    append TextToWrite "\t$FtIncludedInSVtext"
	}
	####### "Breakpoints annotations"
	if {$gcContentText ne ""} {
	    append TextToWrite "\t$gcContentText"
	}
	if {$repeatText ne ""} {
	    append TextToWrite "\t$repeatText"
	}
	####### "Genes-based annotations"
	if {$genesBasedText ne ""} {
	    append TextToWrite "\t$genesBasedText"
	}
	####### "Ranking"
	if {$g_AnnotSV(ranking)} {
	    set rank [SVranking $TextToWrite]
	    if {$g_AnnotSV(rankFiltering) ne ""} {
		if {[lsearch -exact $g_AnnotSV(rankFiltering) $rank] ne -1} {
		    continue
		}
	    }
	    append TextToWrite "\t$rank"
	}

	lappend L_TextToWrite "$TextToWrite"	
    }

    WriteTextInFile [join $L_TextToWrite "\n"] "$outputFile"

    ## Delete temporary file
    file delete -force $tmpFullAndSplitBedFile
    file delete -force $g_AnnotSV(bedFile) ; # => ".formatted.bed" tmp file
    file delete -force $g_AnnotSV(fullAndSplitBedFile)
    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".breakpoints.bed" tmpBreakpointsFile
    file delete -force $tmpBreakpointsFile
    foreach tmpBedFile [glob -nocomplain $g_AnnotSV(outputDir)/*_AnnotSV_inputSVfile.bed] {
	file delete -force $tmpBedFile ; # => a bedfile is present only if "-SVinputFile" is a VCF
    }
    if {[info exists headerFileToRemove]} {
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".header.tsv" BEDinputHeaderFile
	file delete -force $BEDinputHeaderFile
    }

    puts "\n\n...Output columns annotation:"
    regsub -all "\t" $headerOutput "; " t
    puts "\t$t\n"

}

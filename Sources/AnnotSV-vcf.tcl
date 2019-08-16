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



##################################################################
# Prepare 2 new annotation columns for annotated bedfile:
# - nb Hom in the SV
# - nb Htz in the SV
##################################################################


proc DichotomySearch {number OrderedNumberList} {
    # Search for the place of a number into a big ordered number list. 
    # Return the indice after which we can range the number
    # Return -1 if number < [lindex $ListeOrdonnee 0]
    set i 0
    set j [expr {[llength $OrderedNumberList]-1}]
    if {$number<[lindex $OrderedNumberList 0]} {return -1}
    if {$number>=[lindex $OrderedNumberList $j]} {return $j}
    set k [expr {($i+$j)/2}]

    if {$number<[lindex $OrderedNumberList $k]} {set j $k } else {set i $k}
    while {[expr {$j-$i}] > 1} {
        set k [expr {($i+$j)/2}]
	if {$number<[lindex $OrderedNumberList $k]} {set j $k } else {set i $k}
    }
    return $i
}

proc VCFannotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global lPos

    ## VCFs parsing is done only 1 time (After, g_AnnotSV(vcfParsing) is set to "done")
    if {![info exists g_AnnotSV(vcfParsing)]} {
	set g_AnnotSV(vcfParsing) "done"
	# parsing of $g_AnnotSV(vcfFiles): creation of lPos($SVchrom,htz,sample) and lPos($SVchrom,hom,sample)
	puts "\n\n...parsing of VCF file(s) for \"$g_AnnotSV(vcfSamples)\" ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 

	# "eval glob" accept regular expression ("*.vcf) as well as a list of files ("sample1.vcf sample2.vcf.gz"):
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(vcfFiles)] {
	    puts "\t...parsing of $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
	    set iIntersect 0
	    set iSV 0
	    set iLoaded 0
	    set iNotPASS 0
	    set iNotGT 0
	    set chrManipulation 0

	    # Intersect each VCF with the bedfile
	    if {[regexp ".gz$" $vcfF]} {
		regsub ".gz$" [file tail $vcfF] ".tmp.gz" tmpVCFgz
		set tmpVCFgz "$g_AnnotSV(outputDir)/$tmpVCFgz"
		regsub ".gz$" $tmpVCFgz "" tmpVCF

		# VCF should not contain "chr" prefix
		catch {exec gunzip -c $vcfF | grep -c "^chr"} Message
		if {[string index $Message 0] ne 0} {
		    exec gunzip -c $vcfF | sed "s/^chr//" | gzip  > $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf.gz
		    set vcfF $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf.gz
		    set chrManipulation 1
		}

		set f [open "| gzip -cd $vcfF"]	 
		set L_tmpVCF_Header {} ; # Header is used to find "samples" columns names
		while {! [eof $f]} {
		    set L [gets $f]
		    if {[regexp "^#" $L]} {
			set L_tmpVCF_Header $L
		    } else {
			break
		    }
		}
		unset f
		
		WriteTextInFile "$L_tmpVCF_Header" $tmpVCF
		if {[catch {exec gzip $tmpVCF} Message]} {
		    puts "-- VCFannotation --"
		    puts "gzip $tmpVCF"
		    puts $Message
		}
		if {[catch {exec gunzip -c $vcfF | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis} Message]} {
		    if {[regexp -nocase "error" $Message]} {
			WriteTextInFile "$Message" $vcfF.intersect.error
			puts "-- VCFannotation --"
			puts "gunzip -c $vcfF | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis"
			puts "Error: look at $vcfF.intersect.error"
		    }
		}	
	    } else {
		set tmpVCFgz "$g_AnnotSV(outputDir)/[file tail $vcfF].tmp.gz"
		set tmpVCF "$g_AnnotSV(outputDir)/[file tail $vcfF].tmp"
		
		# VCF should not contain "chr" prefix
		catch {exec grep -c "^chr" $vcfF} Message
		if {[string index $Message 0] ne 0} {
		    exec sed "s/^chr//" $vcfF > $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf
		    set vcfF $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf
		    set chrManipulation 1
		}

		set f [open "$vcfF"]	 
		set L_tmpVCF_Header {} ; # Header is used to find "samples" columns names
		while {! [eof $f]} {
		    set L [gets $f]
		    if {[regexp "^#" $L]} {
			set L_tmpVCF_Header $L
		    } else {
			break
		    }
		}
		close $f

		WriteTextInFile "$L_tmpVCF_Header" $tmpVCF
		if {[catch {exec gzip $tmpVCF} Message]} {
		    puts "-- VCFannotation --"
		    puts "gzip $tmpVCF"
		    puts $Message
		}
		if {[catch {exec $g_AnnotSV(bedtools) intersect -a $vcfF -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis} Message]} {
		    if {[regexp -nocase "error" $Message]} {
			WriteTextInFile "$Message" "$g_AnnotSV(outputDir)/[file tail $vcfF].intersect.error"
			puts "-- VCFannotation --"
			puts "$g_AnnotSV(bedtools) intersect -a $vcfF -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis"
			puts "Error: look $g_AnnotSV(outputDir)/[file tail $vcfF].intersect.error"
		    }
		}
	    }

	    if {$chrManipulation} {file delete -force $vcfF}

	    if {[file size $tmpVCFgz.bis] eq 0} {
		# no intersection between the bed and the VCF, no SNV/Indel to load in memory
		puts "\t\t-> no variant in the intersection"
		file delete -force $tmpVCFgz
		file delete -force $tmpVCFgz.bis
		continue
	    } else {
		if {[catch {exec cat $tmpVCFgz.bis | gzip >> $tmpVCFgz} Message]} {
		    puts "-- VCFannotation --"
		    puts "cat $tmpVCFgz.bis | gzip >> $tmpVCFgz"
		    puts $Message
		}
		file delete -force $tmpVCFgz.bis
	    }

	    # load SNV/Indel in memory
	    set f [open "| gzip -cd $tmpVCFgz"]	  
	    set L_samples {}
	    while {![eof $f]} {
		set L [gets $f]

		if {[string range $L 0 5] eq "#CHROM"} {
		    foreach sample $g_AnnotSV(vcfSamples) {
			set i_sample($sample) [lsearch -exact [split $L "\t"] $sample]
			if {$i_sample($sample) ne -1} {lappend L_samples $sample}
		    }
		    continue
		}
		if {[string index $L 0] eq "#" || $L eq ""} {continue}

		incr iIntersect

		# set the variables
		set Ls [split $L "\t"] 
		set chrom  [lindex $Ls 0]
		regsub -nocase "chr" $chrom "" chrom
		set pos    [lindex $Ls 1]

		# Consider only the SNV/indel (not the SV in the VCF file)
		##########################################################
		# Example of SV: 
		# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
		# - Type2: "<INS>", "<DEL>", ...
		# - Type3: complex rearrangements with breakends: "G]17:1584563]"
		set ref [lindex $Ls 3]
		set alt [lindex $Ls 4]
		if {[regexp "<|\\\[|\\\]" $alt]} {incr iSV; continue}; # it is an SV
		set variantLength [expr {[string length $ref]-[string length $alt]}]
		if {[expr {abs($variantLength)}]>$g_AnnotSV(SVminSize)} {incr iSV; continue}; # it is an SV

		set filter [lindex $Ls 6]
		set formatData [lindex $Ls 8]
		set j_GT [lsearch -exact [split $formatData ":"] "GT"]
		if {$j_GT eq -1} {incr iNotGT; continue}
	      
		# keep only variant with FILTER == PASS
		if {$g_AnnotSV(vcfPASS) && $filter ne "PASS"} {incr iNotPASS; continue}
		incr iLoaded
		
		foreach sample $L_samples {
		    set sampleData [lindex $Ls $i_sample($sample)]
		    # set the GT
		    set GTsample [lindex [split $sampleData ":"] $j_GT]
		    set GTsample [split $GTsample "/|\\|"]
		    set GTsampleA [lindex $GTsample 0]
		    set GTsampleB [lindex $GTsample 1]
		    if {$GTsampleA eq "" || $GTsampleB eq ""} {continue}
		    if {$GTsampleA eq "0" && $GTsampleB eq "0"} {continue}
		    if {[lindex $GTsample 0] ne [lindex $GTsample 1]} {set GT "htz"} else {set GT "hom"} 
		    # set the lPos($chrom,$GT,$sample)
		    lappend lPos($chrom,$GT,$sample) $pos
		}
	    }
	    close $f
	    file delete -force $tmpVCFgz
	    if {$iIntersect} {
		puts "\t\t-> $iIntersect variants located in the SV"
	    }
	    if {$iLoaded} {
		puts "\t\t-> $iLoaded variants loaded"
	    }
	    if {$iSV} {
		puts "\t\t-> $iSV SV excluded (considering only SNV/indel from the VCF)"
	    }
	    if {$iNotPASS} {
		puts "\t\t-> $iNotPASS variants excluded beacause of the FILTER value not equal to \"PASS\""
	    }
	    if {$iNotGT} {
		puts "\t\t-> $iNotGT variants excluded because of the absence of GT information"
	    }
	}
    }

    set textToReturn {}
    foreach sample $g_AnnotSV(vcfSamples) {
	# if $chrom not present in the VCF file:
	if {![info exists lPos($SVchrom,htz,$sample)] && ![info exists lPos($SVchrom,hom,$sample)]} {lappend textToReturn "0\t0"; continue}
	# Count for hom variants
	set countHom 0
	if {[info exists lPos($SVchrom,hom,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,hom,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,hom,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,hom,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHom [lsort -unique [lrange $lPos($SVchrom,hom,$sample) $i_first $i_last]]
	    set countHom [llength $L_posHom]
	} 
	# Count for htz variants
	set countHtz 0
	if {[info exists lPos($SVchrom,htz,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,htz,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,htz,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,htz,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHtz [lsort -unique [lrange $lPos($SVchrom,htz,$sample) $i_first $i_last]]
	    set countHtz [llength $L_posHtz]
	}
	lappend textToReturn "$countHom\t$countHtz"
    }

    return "[join $textToReturn "\t"]" 
}


proc VCFsToBED {SV_VCFfiles} {

    global g_AnnotSV
    global VCFheader
    global g_SVLEN

    set SV_BEDfile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d"]_AnnotSV_inputSVfile.bed"
    file delete -force "$SV_BEDfile"
    set VCFheader "" 

    foreach VCFfile $SV_VCFfiles {
	set L_TextToWrite {}
	
	if {[regexp ".gz$" $VCFfile]} {
	    set f [open "| gzip -cd $VCFfile"]	 
	} else {		
	    set f [open "$VCFfile"]
	}	 
	
	set i 0
	while {![eof $f]} {
	    incr i
	    if {$i eq "500000"} {
		WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile
		set L_TextToWrite {}	    
		set i 0
	    }
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    if {[string index $L 0] eq "#" || $L eq ""} {
		if {[regexp "^#CHROM" $L]} {
		    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
		    if {$g_AnnotSV(SVinputInfo)} {
			set VCFheader "SV type\t[join [lrange $Ls 2 end] "\t"]"
		    } else {
			set VCFheader "SV type\t[join [lrange $Ls 3 4] "\t"]\t[join [lrange $Ls 8 end] "\t"]"
		    }
		}
		continue
	    }

	    # Consider only the SV (not the SNV/indel)
	    ##########################################
	    # Example of SV: 
	    # - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG" (length > 50bp)
	    # - Type2: alt="<INS>", "<DEL>", ...
	    # - Type3: complex rearrangements with breakends: alt="G]17:1584563]"
	    set chrom [lindex $Ls 0]
            set pos [lindex $Ls 1]
	    regsub -nocase "chr" $chrom "" chrom
	    set ref [lindex $Ls 3]
	    set alt [lindex $Ls 4]
            if {[regexp "SVLEN=(\[0-9-\]+)" [lindex $Ls 7] match SVLEN]} {set svlen $SVLEN} else {set svlen ""}
	    if {[regexp "END=(\[0-9\]+)" [lindex $Ls 7] match END]} {set end $END} else {set end ""}
	    if {[regexp "SVTYPE=(\[^;\]+)" [lindex $Ls 7] match SVTYPE]} {set svtype $SVTYPE} else {set svtype ""}            
	    if {[regexp "^<" $alt]} {
		# Type2
		if {$end eq ""} {
		    # INS:ME (LINE1, ALU or SVA)
		    set end [expr {$pos+1}]
		} 
		if {[regexp "<TRA>" $alt]} {
		    # Example of a strange VCF format line (translocation):
		    # chr1    63705386   N    .     <TRA>   51      PASS    "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;SVTYPE=TRA;CHR2=chrX;END=478444;SVLEN=0"        GT      ./.   
		    # POS = 63705386 --> sur chrom1
		    # END = 478444   --> sur chromX
		    ## => annotation only of the first breakpoint
		    set end [expr {$pos+1}]
		}
	    } elseif {[regexp "(\\\[|\\\])(\[0-9\XYMT]+):(\[0-9\]+)" $alt match titi chromType3 endType3]} {
		# Type3 		
		if {$chromType3 eq $chrom} {
		    # DEL, DUP, INV...
		    if {$end eq ""} {set end $endType3}
		} else {
		    # TRA (only one breakend wil be annotated with this line)
		    set end [expr {$pos+1}]
		}
		
	    } elseif {[regexp -nocase "^\[ACGTN-\]+$" $ref$alt]} {
		# Type1
		set variantLengthType1 [expr {[string length $alt]-[string length $ref]}]
		if {[expr {abs($variantLengthType1)}]<$g_AnnotSV(SVminSize)} {continue}; # it is an indel
		if {$variantLengthType1>0} {
		    # insertion
		    if {$end eq ""} {set end [expr {$pos+1}]} 
		    if {$svtype eq ""} {set svtype "INS"}
		} else {
		    # deletion
		    if {$end eq ""} {set end [expr $pos-$variantLengthType1]}
		    if {$svtype eq ""} {set svtype "DEL"}
		}
		# Complex SV: AGT>ATTGCATGGACCTGAGTCCCCAAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCGGGGGGGGGG 
		if {![regexp -nocase "^$ref" $alt] && ![regexp -nocase "^$alt" $ref] && $ref ne "-" && $alt ne "-"} {
		    set variantLengthType1 ""
		}
		if {$svlen eq ""} {set svlen $variantLengthType1}
	    } else {
		continue
	    }

	    if {$end <= $pos} {continue}
	    if {$end eq ""} {continue}

	    if {$svtype eq "CNV" || $svtype eq ""} {set svtype $alt}
	
	    
	    if {$g_AnnotSV(SVinputInfo)} {
		lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join [lrange $Ls 2 end] "\t"]"
	    } else {
		lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]"
	    }
	    # Definition of g_SVLEN:
	    # (If not defined here, the variant length can be calculated in AnnotSV-write.tcl for some type of SV)
	    if {$svlen ne ""} {
		set g_SVLEN(${chrom}_${pos}_${end}_$svtype) $svlen
	    }
	}
	close $f
	
	WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile

	if {$VCFheader eq ""} { ; # No header in the VCF input file
	    set length [llength [split "$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]" "\t"]]
	    set VCFheader "SV type\t[join [lrepeat $length ""] "\t"]"
	}

    }

    # Bedfile should be sorted and should not have "chr" in the first column
    # -> This treatment will be done in the 'refGeneAnnotation' proc

    return "$SV_BEDfile"
}

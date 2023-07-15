# This program takes clusters of Ca2+0 ion and determines how many, of these
# ions are in the cluster, how many oxylate counter ions there are how many
# citrate (if applicable) and finally how many water molecules there are and
# the water density within the cluster.

menu main on
mol load psf crys-fin.psf dcd calc.dcd


proc clust { c } {
    set nei_list {}
    set oli_list {}
    set neip {}
#    set nei_cl [atomselect top "same residue as (type CAL or resname OXYL) and within 2.5 of residue $c"]
   set nei_cl [atomselect top "same residue as (type CAL Ni2p or resname OXYL) and within 4.5 of residue $c"]
    set neip [ lsort -unique [$nei_cl get residue]]
#    puts $neip
    set nei_list  [concat $nei_list $neip]
#    puts $nei_list
    set num_nei [ llength $neip ]
    set num_nnei -1

#    puts $num_nei
#    puts $num_nnei
    while {$num_nnei != 0 } {
        set new_nei {}
        set clnei "$nei_list"
#	puts $clnei
#        set nnei_cl [atomselect top "((same residue as (type CAL or resname OXYL)) and (within 2.5 of residue $clnei)) and (not residue $clnei)"]
        set nnei_cl [atomselect top "((same residue as (type CAL Ni2p or resname OXYL)) and (within 4.5 of residue $clnei)) and (not residue $clnei)"]


# Check cluster
        set new_nei [ lsort -unique [$nnei_cl get residue]]
#	puts $new_nei

# set criterion for while loop
        set num_nnei [ llength $new_nei ]

        set nei_list [concat $nei_list $new_nei ]
        set oli_list [ lsort -unique $nei_list ]

        set num_nei [ llength $oli_list]

     
#        puts $num_nei
#        puts $num_nnei
        $nnei_cl delete

        array unset $new_nei
    }

#    puts $nei_list

    set oli_list [ lsort -unique $nei_list ]
    $nei_cl delete
#    puts $oli_list

    return $oli_list
}

proc dimensions { s } {
    set sel [atomselect top "residue $s"]
    set pts [measure minmax $sel]
    set finalcoor {}

    foreach { c1 c2 } $pts { break }
    set c11 $c1
    set c22 $c2

    foreach { x0 y0 z0 } $c11 { break }
    set cx0 $x0
    set cy0 $y0
    set cz0 $z0

    foreach { x1 y1 z1 } $c22 { break }
    set cx1 $x1
    set cy1 $y1
    set cz1 $z1

    $sel delete
    set finalcoor "$cx0 $cy0 $cz0 $cx1 $cy1 $cz1"
    return $finalcoor
}

proc trans { m1 m2 cond val } {
    set m [concat $m1 $m2]
    set clu [atomselect top "residue $m"]

    set vec {0 0 0}
    set vec [lreplace $vec $cond $cond $val]

    set d0 [lindex $vec 0]
    set d1 [lindex $vec 1]
    set d2 [lindex $vec 2]
    set disp "$d0 $d1 $d2"

    $clu moveby "$disp"
    $clu delete
    return
}


set Caf [open Ca-cluster.xvg w]
set Oxf [open Ox-cluster.xvg w]
set Citf [open Cit-cluster.xvg w]
set Watf [open Wat-cluster.xvg w]

puts $Caf "time Cluster Ca2+ Wat:Ca2+"
puts $Oxf "time Cluster Oxyl Wat:Oxyl"
puts $Citf "time Cluster Citr Wat:Citr"
puts $Watf "time Cluster Wat"

#set Cluster [atomselect top "type CAL or resname OXYL"]
set Cluster [atomselect top "type CAL Ni2p or resname OXYL"]
set Calist [lsort -unique [$Cluster get residue]]
set nf [molinfo top get numframes]

## determine sizes in case cluster is split along PBC
#set all [atomselect top all]
#set a [ $all get residue]
#$all delete
#array unset a
#set box [ dimensions $a ]
#foreach { xa ya za xb yb zb xl yl zl } $box { break }
#set xbmin $xa
#set ybmin $ya
#set zbmin $za
#set xbmax $xb
#set ybmax $yb
#set zbmax $zb
#set xbox $xl
#set ybox $yl
#set zbox $zl

set xbmin -100.61
set ybmin -100.83
set zbmin -100.76
set xbmax 101.28
set ybmax 101.37
set zbmax 99.97
set xbox 201.89
set ybox 202.20
set zbox 200.74
for {set i 0} {$i < $nf} {incr i} {
    molinfo top set frame $i
    set Cadone {}
    foreach d $Calist {
	set pres [lsearch -all -inline $Cadone $d]
	set lpres [ llength $pres]

#        puts $pres
#        puts $lpres

	if { $pres == {} } {

# determine clusters first time
	    set lists [ clust $d ]
#	    puts $lists

	    set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
	    set neil_Wat [lsort -unique [$wat_cl get residue]]	        
#	    puts $neil_Wat
    
#determine whether clusters are split along PB, if so translate them
	    set cludim [ dimensions $lists ]
#            puts $cludim
	    foreach { xa ya za xb yb zb } $cludim { break }
	    set xcmin $xa
	    set ycmin $ya
	    set zcmin $za
	    set xcmax $xb
	    set ycmax $yb
	    set zcmax $zb


	    set dise1 [ expr { abs($xcmin - $xbmin) } ]
	    set dise2 [ expr { abs($xbmax - $xcmax) } ]

#	    puts $dise1
#	    puts $dise2
	    if { $dise1 < 4.5 } {
	        trans $lists $neil_Wat 0 $xbox
	        set lists [ clust $d ]
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]          
                set neil_Wat [lsort -unique [$wat_cl get residue]]
	    }

	    if { $dise2 < 4.5 } {
	        trans $lists $neil_Wat 0 [expr {-$xbox}]
                set lists [ clust $d ]
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
                set neil_Wat [lsort -unique [$wat_cl get residue]]
	    }
	    
	    set dise1 [ expr { abs($ycmin - $ybmin) } ]
            set dise2 [ expr { abs($ybmax - $ycmax) } ]

#            puts $dise1
#            puts $dise2
            if { $dise1 < 4.5 } {
                trans $lists $neil_Wat 1 $ybox
                set lists [ clust $d ]
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
                set neil_Wat [lsort -unique [$wat_cl get residue]]
            }
            if { $dise2 < 4.5 } {
                trans $lists $neil_Wat 1 [expr {-$ybox}]
                set lists [ clust $d ]
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
                set neil_Wat [lsort -unique [$wat_cl get residue]]
            }

            set dise1 [ expr { abs($zcmin - $zbmin) } ]
            set dise2 [ expr { abs($zbmax - $zcmax) } ]

#            puts $dise1
#            puts $dise2
            if { $dise1 < 4.5 } {
                trans $lists $neil_Wat 2 $zbox
                set lists [ clust $d ]
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
                set neil_Wat [lsort -unique [$wat_cl get residue]]
            }
            if { $dise2 < 4.5 } { 
                trans $lists $neil_Wat 2 [ expr {-$zbox}]
                set lists [ clust $d ]
# 		puts $lists
                set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]
                set neil_Wat [lsort -unique [$wat_cl get residue]]
#		puts $neil_Wat
            }

	    set Cadone [ concat $Cadone $lists ]

#	    puts $Cadone

# Now that whole cluster has been determined we will count number in each
# species and the densities

	    set neil_Ca [atomselect top "type CAL and residue $lists"]
	    set neil_Ox [atomselect top "resname OXYL and residue $lists"]
	    set neil_Cit [atomselect top "type Ni2p and residue $lists"]
	    set wat_cl [atomselect top "same residue as water and within 4.5 of residue $lists"]

#	    puts [lsort -unique [$neil_Ca get residue]]
#	    puts [lsort -unique [$neil_Ox get residue]]


	    set numCa [llength [lsort -unique [$neil_Ca get residue]]]
	    set numOx [llength [lsort -unique [$neil_Ox get residue]]]
	    set numCit [llength [lsort -unique [$neil_Cit get residue]]]
	    set numWat [llength [lsort -unique [$wat_cl get residue]]]

	    set numCa [ expr {$numCa+0.0} ]
            set numOx [ expr {$numOx+0.0} ]
            set numCit [expr {$numCit+0.0} ]
            set numWat [expr {$numWat+0.0} ]

# Determine ratio number of water molecules around Ca2+ cluste
	    if {
               $numCa != 0
               && [expr { $numOx + $numCit}] > 0 
       } then {
            set Caratiwat [ expr { $numWat/$numCa } ]
            set Oxratiwat [ expr { $numWat/$numOx } ]
	    set Citratiwat [ expr { $numWat/$numCit } ]
	    # Now output data
	    puts $Caf "$i $numCa $Caratiwat"
	    puts $Oxf "$i $numOx $Oxratiwat"
	    puts $Citf "$i $numCit $Citratiwat" 
	    puts $Watf "$i $numWat"
    }

            $neil_Ca delete
            $neil_Ox delete
            $neil_Cit delete
            $wat_cl delete
	    array unset $lists
        } 
    }
}

$Cluster delete
#$all delete
exit

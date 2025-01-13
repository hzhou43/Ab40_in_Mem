#load trajectory
mol load psf ../step5_input.psf pdb ../step5_input.pdb
mol addfile ../step7_production.dcd first 0 last -1 step 5 waitfor all

#output file
set outfile [open tilt_angle.txt w]

#getting the segnames of all beta strands
set pro [atomselect top "protein"]
set seg [lsort -unique [$pro get segname]]

#starting residue of the beta strand
set start 16

#ending residue of the beta strand
set end 21

#total number of frames
set nf [molinfo top get numframes]

for {set j 0 } { $j <  $nf} {incr j} {
	for {set i 0 } { $i <  3} {incr i} {
		#selecting C and N atoms 
		set sel1 [atomselect top "segname [lindex $seg $i] and resid $start and name C " frame $j]
		set sel2 [atomselect top "segname [lindex $seg $i] and resid [expr {$start+1}] and name N " frame $j]
		set sel3 [atomselect top "segname [lindex $seg $i] and resid [expr {$end-1}] and name C " frame $j]
		set sel4 [atomselect top "segname [lindex $seg $i] and resid $end and name N" frame $j]

		set sel5 [atomselect top "segname [lindex $seg $i+1] and resid $start and name C " frame $j]
                set sel6 [atomselect top "segname [lindex $seg $i+1] and resid [expr {$start+1}] and name N " frame $j]
                set sel7 [atomselect top "segname [lindex $seg $i+1] and resid [expr {$end-1}] and name C" frame $j]
                set sel8 [atomselect top "segname [lindex $seg $i+1] and resid $end and name N" frame $j]


		#getting xyz coordinates of C and N atoms
		set cord1 [lindex [$sel1 get {x y z}] 0]
		set cord2 [lindex [$sel2 get {x y z}] 0]
		set cord3 [lindex [$sel3 get {x y z}] 0]
		set cord4 [lindex [$sel4 get {x y z}] 0]
		set cord5 [lindex [$sel5 get {x y z}] 0]
                set cord6 [lindex [$sel6 get {x y z}] 0]
                set cord7 [lindex [$sel7 get {x y z}] 0]
                set cord8 [lindex [$sel8 get {x y z}] 0]

		#midpoint of C-N bonds of two consecutive strands
                set A [vecscale 0.5 [vecadd $cord1 $cord2]]
		set B [vecscale 0.5 [vecadd $cord3 $cord4]]
		set C [vecscale 0.5 [vecadd $cord5 $cord6]]
		set D [vecscale 0.5 [vecadd $cord7 $cord8]]

		#calculation of the angle between the two vectors
		set AB [vecsub $B $A]
		set CD [vecsub $D $C]

		set dot_product [vecdot $AB $CD]


		set mag_AB [expr {sqrt([vecdot $AB $AB])}]
		set mag_CD [expr {sqrt([vecdot $CD $CD])}]
		set cos_theta [expr {$dot_product / ($mag_AB*$mag_CD)}]

		set angle_radians [expr {acos($cos_theta)}]
		set angle_degrees [expr {$angle_radians * 180.0 / acos(-1)}]
	
		puts $outfile "$angle_degrees"
	}
}

close $outfile

exit

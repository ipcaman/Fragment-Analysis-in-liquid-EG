# tcl script to provide input data for SOM
# Steps - 1: To build different fragments of EG molecules in liquid state; trimers, tetramers and pentamers 
#         2: For the corresponding fragment save its dihedral angle values, frame information and indexes of molecules involved
#         3: Save the coordinates of all the atoms forming the selected fragment 

# procedure to dump frame and molecular indexes information into one file, dihedral angles and coordinates of all the atoms forming the fragment in separate files
# these files (or say data) will directly serve as MATLAB input file
# input required for the following procedure: indexes of the molecules forming the fragment along with the respective frame.  
proc dihed {frame args} {
    set mer {}
    # for trimers 
    if {[llength [lindex $args 0]] == 3} {
    set fp1 [open "trimers_frames_mol.dat" "a"]
    set fp2 [open "trimers_dihedrals.dat" "a"]
    set fp3 [open "trimers_pax.dat" "a"]
    #for tetramers
    } elseif {[llength [lindex $args 0]] == 4} {
    set fp1 [open "tetramers_frames_mol.dat" "a"]
    set fp2 [open "tetramers_dihedrals.dat" "a"]
    set fp3 [open "tetramers_pax.dat" "a"]
    #for tetramers
    } elseif {[llength [lindex $args 0]] == 5} {
    set fp1 [open "pentamers_frames_mol.dat" "a"]
    set fp2 [open "pentamers_dihedrals.dat" "a"]
    set fp3 [open "pentamers_pax.dat" "a"]
    #for pentamers
    } elseif {[llength [lindex $args 0]] == 6} {
    set fp1 [open "hexamers_frames_mol.dat" "a"]
    set fp2 [open "hexamers_dihedrals.dat" "a"]
    set fp3 [open "hexamers_pax.dat" "a"]
    }

    foreach oxy $args {
        
        foreach  oxy_index $oxy {
            set atom1 {}
            set o {}
            set h {}
            set c {}
            set mol_ [atomselect top "same fragment as index $oxy_index"] 
            set mol [$mol_ get index]
            $mol_ delete
            foreach q $mol {
                set dummy_check_ [atomselect top "index $q"] 
                set dummy_check [$dummy_check_ get type]
                $dummy_check_ delete
                if {$dummy_check == "O"} {
                    lappend o $q
                } elseif {$dummy_check == "H"} {
                    lappend h $q
                } else {
                    lappend c $q
                }
            }


            lappend atom1 [lindex $o 0]
            lappend atom1 [lindex $c 1]
            lappend atom1 [lindex $c 0]
            lappend atom1 [lindex $o 1]
            
            set value1 [measure dihed $atom1]
            set OCCO [format "%.2f" $value1] 
            lappend mer [format "%.2f" [expr {abs($OCCO)}] ]  
        }
    }

    #to get the atomic coordinates
    set pax [[atomselect top "same fragment as index [lrange [lindex $args 0] 0 end]"] get {x y z}]   
    #dumping frame and molecular indexes in file fp1 
    puts $fp1 "$frame $args"
    #fp2 contains dihedral angle information for the corresponding fragment
    puts $fp2 "$mer"

    foreach pac $pax { 
       #dumping atomic coordinates in fp3 
       puts $fp3 "[format "%.2f" [lindex $pac 0]] [format "%.2f" [lindex $pac 1]] [format "%.2f" [lindex $pac 2]]" 
    }    
    close $fp1
    close $fp2
    close $fp3
    return 
}

#above procedure ends here



#selecting the molecules within 3 Angstrom on each side of the box
set index_ot_ [atomselect top "type O and same fragment as (x > -6 and x < 25 and  y > -6 and y < 26  and z  > -4 and z < 27)"]
set index_o [$index_ot_ get index]
$index_ot_ delete

set element {}
#element list contains the indexes of oxygen atoms for all 100 EG molecules within the box
for {set x 0} {$x < [llength $index_o]} {incr x 2} {
    lappend element [lindex $index_o $x]
}

# to get total number of frames
set nf [molinfo top get numframes]

for {set f 0} {$f < $nf} {incr f 5} {
    molinfo top set frame $f
    puts "I am working on frame $f; Have patience and go for coffee"
    set   h_bonded {}
    set   h_bonded_act {}
    foreach {indoa1 indoa2} $index_o {
        set h_bonded_dimer {}
        set h_bonded_dimer_act {}
        if {($indoa1  > 13399) || ($indoa1  < 13200)} {
            set n [expr {$indoa1 % 1000}]
            set indoa1_orig [expr {$n+13000}]
        } else {
            set indoa1_orig $indoa1
        }

        lappend h_bonded_dimer $indoa1_orig
        lappend h_bonded_dimer_act $indoa1

        foreach {indob1 indob2} $index_o {   
            if {($indob1  > 13399) || ($indob1  < 13200)} {
                set n [expr {$indob1 % 1000}]
                set indob1_orig [expr {$n+13000}]
            } else {
                set indob1_orig $indob1
            } 

            if {$indob1_orig != $indoa1_orig} {
                set h_indb1b2_ [atomselect top "type H and within 1.2 of index $indob1 $indob2"]
                set h_indb1b2 [$h_indb1b2_ get index]
                $h_indb1b2_ delete
    
                set h_inda1a2_ [atomselect top "type H and within 1.2 of index $indoa1 $indoa2"]
                set h_inda1a2 [$h_inda1a2_ get index]
                $h_inda1a2_ delete
                   
                #To check the hydrogen bonds
                foreach oa "$indoa1 $indoa2" ha $h_inda1a2 {
                    foreach hb $h_indb1b2 ob "$indob1 $indob2" {
                      
                        # from a to b
                        set bond_distab {}
                        set bond_angleab {}
                        lappend bond_distab $oa
                        lappend bond_distab $hb
                        lappend bond_angleab $oa
                        lappend bond_angleab $hb
                        lappend bond_angleab $ob
  
                        # from b to a
                        set bond_distba {}
                        set bond_angleba {}
                        lappend bond_distba $ob
                        lappend bond_distba $ha
                        lappend bond_angleba $ob
                        lappend bond_angleba $ha
                        lappend bond_angleba $oa
  
  
                        set rab [measure bond $bond_distab]
                        set rba [measure bond $bond_distba]
                        set ang_abb [measure angle $bond_angleab]
                        set ang_bab [measure angle $bond_angleba]
                                 
        
                        if {(($ang_abb >= 140) && ($rab <= 2.40)) || (($ang_bab >= 140) && ($rba <= 2.40))} {
                            lappend h_bonded_dimer $indob1_orig
                            lappend h_bonded_dimer_act $indob1
                            break
                        }
                    }
                    if {(($ang_abb >= 140) && ($rab <= 2.40)) || (($ang_bab >= 140) && ($rba <= 2.40))} {
                        break
                    }
  
                }   
              
            }
        }
        lappend h_bonded $h_bonded_dimer
        lappend h_bonded_act $h_bonded_dimer_act
    }



#procedure to map the indexes of atoms outside the box to their cooresponding indexes within the box; required to avoid double-counting
  proc act_orig {act_ind} {
    global ind_orig
    set orig_ind {}
    foreach ind $act_ind {
        if {($ind   > 13399) || ($ind   < 13200)} {
        set nn [expr {$ind  % 1000}]
        set ind_orig [expr {$nn+13000}]
        } else {
            set ind_orig $ind 
        }
        
    }
    return $ind_orig
  }

  global ind_orig




# Picking up the Trimer fragment
  set trimer {}
  set trimer_act {}


  foreach cluster $h_bonded cluster_act $h_bonded_act {
        set m [llength $cluster]
        if {[llength $cluster] > 2} {
            for {set i 1} {$i < $m} {incr i} {
                for {set j [expr {$i+1}]} {$j < $m} {incr j} {
                    lappend trimer  [lsort "[lindex $cluster 0]  [lindex $cluster $i]   [lindex $cluster $j]"]
                    lappend trimer_act  "[lindex $cluster_act 0]  [lindex $cluster_act $i]  [lindex $cluster_act $j]"
                }
            }
        }
  }
  
  set uniq_trimer [lsort -unique $trimer] 
  set trimer_ans {}

# Taking only the unique ones
  foreach tri $uniq_trimer {
    set ind_trimer [lsearch -exact $trimer $tri]
    lappend trimer_ans [lindex $trimer_act [lindex $ind_trimer 0]]
  }
  puts "trimers = [llength $trimer_ans]"

# calling the procedure dihed
  foreach ans1 $trimer_ans {
     dihed $f  "$ans1"
  }



# Building Tetramers from trimers by adding of one EG molecule in all possible ways 
  set tetramer {}
  set tetramer_act {}            
  foreach tri1 $trimer_ans {
    set a_act [lindex $tri1 0]
    set a [act_orig $a_act]
    set b_act [lindex $tri1 1]
    set b [act_orig $b_act]
    set c_act [lindex $tri1 2]
    set c [act_orig $c_act]
    
#tetramers from a 
    set apair [lindex $h_bonded [lsearch -exact $element $a]]
    set apair_act [lindex $h_bonded_act [lsearch -exact $element $a_act]]

    set lena [llength $apair]
    if {$lena > 2} {
        if {[llength $apair] == [llength $apair_act]} {
            foreach k $apair k_act $apair_act {
                if {($k != $a) && ($k != $b) && ($k != $c)} {
                    lappend tetramer [lsort "$a $b $c $k"]
                    lappend tetramer_act "$a_act $b_act $c_act $k_act"
                }
            }
        }
    }

                         
#tetramer from b
    set bpair [lindex $h_bonded [lsearch -exact $element $b]]
    set bpair_act [lindex $h_bonded_act [lsearch -exact $element $b_act]]
    set lenb [llength $bpair]
    if {$lenb > 2} {
        if {[llength $bpair] == [llength $bpair_act]} {
            foreach l $bpair l_act $bpair_act {
                if {($l != $a) && ($l != $b) && ($l != $c)} {
                    lappend tetramer [lsort "$a $b $c $l"]
                    lappend tetramer_act "$a_act $b_act $c_act $l_act"
                }
            }
        }
    }

#tetramers from c

    set cpair [lindex $h_bonded [lsearch -exact $element $c]]
    set cpair_act [lindex $h_bonded_act [lsearch -exact $element $c_act]]
    set lenc [llength $cpair]

    if {$lenc > 2} {
        if {[llength $cpair] == [llength $cpair_act]} {
            foreach o $cpair o_act $cpair_act {
                if {($o != $a) && ($o != $b) && ($o != $c)} {
                    lappend tetramer [lsort "$a $b $c $o"]
                    lappend tetramer_act "$a_act $b_act $c_act $o_act"
                }
            }
        }
    }
               
  }
                            
#deleting the repeating ones                         
  set uniq_tetra [lsort -unique $tetramer] 
  set tetra_ans {}

  foreach tet $uniq_tetra {
    set ind_tet [lsearch -exact $tetramer $tet]
    lappend tetra_ans [lindex $tetramer_act [lindex $ind_tet 0]]
  }
  puts "tetramers = [llength $tetra_ans]"

# calling the procedure dihed
  foreach ans2 $tetra_ans {
    dihed $f  "$ans2"
   }



# Building Pentamers from tetramers by adding of one EG molecule in all possible ways 
  set pentamer {}
  set pentamer_act {}

  foreach tetra1 $tetra_ans {
    set a1_act [lindex $tetra1 0]
    set a1 [act_orig $a1_act]
    set b1_act [lindex $tetra1 1]
    set b1 [act_orig $b1_act]
    set c1_act [lindex $tetra1 2]
    set c1 [act_orig $c1_act]
    set d1_act [lindex $tetra1 3]
    set d1 [act_orig $d1_act]
                    
    #puts "petamers from a" 
    set a1pair [lindex $h_bonded [lsearch -exact $element $a1]]
    set a1pair_act [lindex $h_bonded_act [lsearch -exact $element $a1_act]]
    set lena1 [llength $a1pair]
    if {$lena1 > 2} {
        if {[llength $a1pair] == [llength $a1pair_act]} {
            foreach q $a1pair q_act $a1pair_act {
                if {($q != $a1) && ($q != $b1) && ($q != $c1) && ($q != $d1)} {
                    lappend pentamer_act "$a1_act $b1_act $c1_act $d1_act $q_act"
                    lappend pentamer [lsort "$a1 $b1 $c1 $d1 $q"]
                }
            }           
        }
    }
                    
#puts "pentamers from b"
    set b1pair [lindex $h_bonded [lsearch -exact $element $b1]]
    set b1pair_act [lindex $h_bonded_act [lsearch -exact $element $b1_act]]
    set lenb1 [llength $b1pair]
    if {$lenb1 > 2} {
        if {[llength $b1pair] == [llength $b1pair_act]} {
            foreach l1 $b1pair l1_act $b1pair_act {
                if {($l1 != $a1) && ($l1 != $b1) && ($l1 != $c1) && ($l1 != $d1)} {
                    lappend pentamer_act "$a1_act $b1_act $c1_act $d1_act $l1_act"
                    lappend pentamer [lsort "$a1 $b1 $c1 $d1 $l1"]
                }
            }           
        }
    }

    #puts "pentamers from c"
    set c1pair [lindex $h_bonded [lsearch -exact $element $c1]]
    set c1pair_act [lindex $h_bonded_act [lsearch -exact $element $c1_act]]
    set lenc1 [llength $c1pair]
    if {$lenc1 > 2} {
        if {[llength $c1pair] == [llength $c1pair_act]} {
            foreach k1 $c1pair k1_act $c1pair_act {
                if {($k1 != $a1) && ($k1 != $b1) && ($k1 != $c1) && ($k1 != $d1)} {
                    lappend pentamer_act "$a1_act $b1_act $c1_act $d1_act $k1_act"
                    lappend pentamer [lsort "$a1 $b1 $c1 $d1 $k1"]
                }
            }           
        }
    }

    #puts "pentamers from d"
    set d1pair [lindex $h_bonded [lsearch -exact $element $d1]]
    set d1pair_act [lindex $h_bonded_act [lsearch -exact $element $d1_act]]
    set lend1 [llength $d1pair]
    if {$lend1 > 2} {
        if {[llength $d1pair] == [llength $d1pair_act]} {
            foreach m1 $d1pair m1_act $d1pair_act {
                if {($m1 != $a1) && ($m1 != $b1) && ($m1 != $c1) && ($m1 != $d1)} {
                    lappend pentamer_act "$a1_act $b1_act $c1_act $d1_act $m1_act"
                    lappend pentamer [lsort "$a1 $b1 $c1 $d1 $m1"]
                }
            }           
        }
    }
  }

#taking only the uniq ones
  set uniq_penta [lsort -unique $pentamer] 
  set penta_ans {}

  foreach pent $uniq_penta {
    set ind_pent [lsearch -exact $pentamer $pent]
    lappend penta_ans [lindex $pentamer_act [lindex $ind_pent 0]]
  }

  puts "pentamers = [llength $penta_ans]"

  # calling the procedure dihed
  foreach ans3 $penta_ans {
    dihed $f  "$ans3"
  }
}
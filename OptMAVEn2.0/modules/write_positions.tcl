source /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/modules/VMD_FUNCTIONS.tcl
#set posFile [open "positions.dat" "a"]
set pos [getAntigenPosition [atomselect 0 "all"] [atomselect 0 "resid 306 307 342 343 344 350 353 355 388 391 392 393 395"] 0]
puts $pos
#close $posFile
exit


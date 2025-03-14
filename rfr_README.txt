//This is a document explaining the code contributed by Rafael Ruiz(rfruiz)
----------------------------------------------------------------------------------


######## DBread.h
----------------------------------------------------------------------------------

### APV key, value map creation
	readDBfile() now also fills the keys(module Number, axis, and position along the axis),value(vtpcrate, fiber, and adc_ch) apvInfoMap that is ultimately used to obtain the defining APV values corresponding with ECalBins. 


### Static APV chanMap data is now in db_FT_local.dat
	The keys and values are read from the new lines that I pasted into db_FT_local.dat.

### How the chanMap data is obtained
	The text was copied from chanMapAPVinfo.txt, the file from GetchanMap.cxx, which parses a desired db file as reference(in this case used db_sbs.gemFT.dat).

	So if chanMap data was to change, the script could be rerun and the new data could replace the current text in db_FT_local.dat




####### generate_ROI_APVlist.C
----------------------------------------------------------------------------------
## New lines demarcated by //RREDIT comments which I can remove later
	I just want the additions to be easier to notice/follow/understand


### new container ECalBinAPVinfo
	I made a new container ECalBinAPVinfo that stores, for each ECalBin, the relevant APV info values in a set.

### filling ECalBinAPVinfo
	ECalBinAPVinfo is filled within the same loops that fill the U and V strip counts
	by 2 for loops, one for U strips and the other for V strips

	Each loop already has the ECalBin and modID from outer loops and the axis is either u or v
	so the only remaining key is "pos" which I just get from 
	the integer division: stripID/128 
	which rounds to the corresponding APV since each should have 128 strips.

	at the end of each of these loops(overt the strips in the u and v axes), if the current key combination exists in the apvInfoMap(if they are not found they are added to a set of missingkeys),
	the APV info values are obtained with the keys and are placed in the set of the current ECalBin APVs

### output file(EcalBinToAPVmap.txt)
writes out the APV values for all ECalBins




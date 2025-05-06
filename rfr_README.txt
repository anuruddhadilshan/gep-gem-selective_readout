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




####### showGEMstripHit_for_ecalbin.cpp
----------------------------------------------------------------------------------
# README for GEM Strip Visualization Tool



##### TLDR

#1)Run with .x showGEMstripHit_for_ecalbin.cpp

You’ll be prompted to:

#2) Enter one or more ECal bin numbers (or press return for all).
	#press q to quit
	# choose one or more bins then return to makes canvas for GEM ROIs for ECalbins
		#ex. "0 90 183" makes canvas for GEM ROIs for ECalbins 0 90 and 183
	#press return to plot for all ECal bins

#3) Choose whether to save the output as a PDF.(y/n)
	



##Details
--------------------------

## Overview
This ROOT-based C++ script visualizes active U/V GEM detector strips for specified ECal bins using Region of Interest (ROI) mappings. It reads static database and ROI configuration files, maps projected hit regions to GEM strip numbers, and draws overlayed strip lines for each layer/module on a multi-panel ROOT canvas.

## Key Features
- Parses detector database and ROI files.
- Maps physical ROIs to U/V strip numbers per GEM module.
- Visualizes active strips per GEM layer for one or multiple ECal bins.
- Outputs to interactive ROOT canvases or saves to PDFs.
- Annotates strip geometry, ROI rectangles, and module boundaries.

## Dependencies
- ROOT libraries (`TCanvas`, `TLine`, `TH2I`, `TLatex`, `TMarker`, `TRandom`, etc.)
- Local project files:
  - `DBread.h/cxx`: Loads GEM geometry and configuration.
  - `ROIread.h/cxx`: Parses ECal ROI definitions.
  - `GEMModROItoStrips.h`: Converts ROI coordinates to physical strip sets.
  - `GetGemInfoMap.cxx`: Returns a map of GEM module info (positions, angles, sizes, etc.).

## Usage
### 1. Compilation
If compiled:
```bash
g++ -o showgem showgemhit_for_ecalbin.C $(root-config --cflags --libs)
```

### 2. Execution (with ROOT CINT or compiled)
```bash
root -l -b -q 'showgemhit_for_ecalbin.C("db_FT_local.dat", "ROI_GEP3_FT_1.txt")'
```

### 3. Interactive Prompts
You’ll be prompted to:
- Enter one or more ECal bin numbers (or press return for all).
- Choose whether to save the output as a PDF.

## Output
- For each ECal bin: an 8-panel canvas with each GEM layer.
- Each panel shows:
  - Dummy histogram with module axes.
  - Strip lines (U = cyan, V = pink).
  - Green bounding box for ROI.
  - Module boundary rectangle.
- PDF saved as `GEM_ECalBin<bin>.pdf` if selected.

## Coordinate Systems
- GEM internal coordinates: +X down, +Y left.
- Strip angles are measured w.r.t. +X.
- Plots are rotated to match global canvas orientation (X → Y, Y → X).

## Customization Notes
- `drawGEMStrip()` adjusts spacing and visibility (e.g., draw every 20th strip).
- Offset corrections applied via `getOffset(mod, axis)` for module alignment.
- Canvas ID and subpad layout configurable (currently 4x2).

## Troubleshooting
- **Web GUI errors**: ROOT may attempt browser-based rendering; disable or switch to file output.
- **Null canvases**: Likely due to missing ROI data for that ECal bin.
- **Strip overlap or misalignment**: Check offset constants and coordinate flipping logic.

## Author
Rafael Ruiz

## License
This script is provided as part of the GEp-V visualization toolkit at Jefferson Lab.

---
Let me know if you'd like to auto-generate `.dox`/Doxygen-compatible comments for each function.


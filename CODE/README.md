cleanupExcel.sh = run this before anything else if you have edited xlsx files - deletes temporary files starting with "~"
readSnowpitXLS.py = reads in all snowpits, and currently just produces a summary of location, depth, density, SWE
readSnowpitXLS_script.py = script version for single file, used for debugging

To rerun after editing xlsx files:
>>./cleanupExcel.sh
>>python readSnowpitXLS.py

NOTES

Pit sheets were entered into xlsx files in the field.
Removed asteriks on observations that made data import difficult.  These notes are now in Comments section.
Added missing coordinates to B01 and B04, from list of originally planned locations (Downselect_50_pts_GM_Mar2019.txt)


  TO_DO
Add date to summary file
extend readSnowpitXLS to write csv files for all variables in SnowEx17 format

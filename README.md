#SAMCompare

The script `compare.py` creates a dictionary of SAM headers objects that hold all of the records for each header. Each record contains a few significant fields: FLAG, MAPQ, RNAME, POS, AS, and NM. By comparing the header objects, SAMCompare helps to answer the following questions:
	
- Alignment Differences
	- Lists FLAG, MAPQ, RNAME, POS, AS, NM for each alignment
	- Statistics on differences in alignments
		- Is the mapped position different?
		- Has the number if alignments per header changed?
- Statitics for each SAM file
	- Headers 
	- Alignments
	- Alignments per Header
	- Records per Header
	- Records
	- Unmapped Alignments

##Version

Version 1.0 supports only paired end reads. Single end reads inclusion will come with later versions.

##Setup
- Make sure latest version of python is installed (Python 3.5.1)
- `$ python3               # check version with terminal command`
- Clone or download repository
- See Manual for running instructions and parameters

##Credits

Written by Stephen Petrides

The MIT License

Copyright (C) 2016 Stephen Petrides

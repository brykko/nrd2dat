# README #

nrd2dat is a utility for converting Neuralynx raw data ".nrd" files to the flat binary ".dat" format used by spike sorting packages such as Klusta, Kilosort and Phy.  It uses buffered IO streams for fast reading and writing of large data files.  The speed of your hard drive will likely be the limiting factor for performance.  A 300 GB  .nrd file is typically processed in ~1 hour on my (fairly slow) system.  Error-checking is performed to ensure records have been read correctly from the raw data file.

### USAGE ###

nrd2dat.exe "SourceNlxRawDataFile.nrd" "ChannelMap.txt"

The Neuralynx .nrd files include all available channels, given the system's input board configuration.  This converter requires a channel map indicating the numbers of the AD channels to be extracted.  Samples from these channels will be written to the .dat file in the order specified in the channel map.  The channel map must be a text file with the number of a single AD channel on each line. Any lines that are empty or  prefixed with '%' will be ignored.  See below for an example channel map.

The Neuralynx raw samples are stored with 24-bit precision, with much larger dynamic range (+- 132 mV) than is typically required for representing the frequency band containing spikes (~600 6000 Hz).  To reduce the loss of precision in converting from 24-bit to 16-bit, the signal is high-pass-filtered at a cutoff of 100 Hz to remove high-amplitude, low-frequency signal components.

### OUTPUT FILES ###
The output is a pair of .dat files whose names begin with the base name of the input raw data file.  Both files are headerless.  The first file <basename>_samples.dat contains a stream of 16-bit signed integers ordered first by AD channel (as specified in the channel map) and second by sample.  The second file <basename>_timestamps.dat contains a stream of 64-bit unsigned integers storing the timestamp of each record.


###EXAMPLE CHANNEL MAP (corresponding to tetrodes 1, 2 and 3 on the HS54 headstage):###

% TT1  
40  
41  
42  
43  

%TT2  
55  
54  
53  
52  

%TT3  
15  
14  
13  
12  

*/


### How do I get set up? ###

Clone the repository and open the solution file "nrd2dat2.sln" in Microsoft Visual Studio.  Build the solution to generate the executable.
# README #

nrd2dat is a utility for converting Neuralynx raw data ".nrd" files to the flat binary ".dat" format required used in packages such as Klusta, Kilosort and Phy.  It uses buffered IO streams for fast reading and writing of large data files.

The Neuralynx .nrd files include all available channels, given the system's input board configuration.  This converter requires a channel map indicating the numbers of the AD channels to be extracted.  Samples from these channels will be written to the .dat file in the order specified in the channel map.  The channel map must be a text file with the number of a single AD channel on each line. Any lines that are empty or  prefixed with '%' will be ignored.  See below for an example channel map.

The output is a pair of .dat files whose names begin with the base name of the input raw data file.  Both files are headerless.

The first file <basename>_samples.dat contains a stream of 16-bit signed integers ordered first by AD channel (as specified in the channel map) and second by sample.

The second file <basename>_timestamps.dat contains a stream of 64-bit unsigned integers storing the timestamp of each record.


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

Usage:
nrd2dat.exe "SourceNlxRawDataFile.nrd" "ChannelMap.txt"

### How do I get set up? ###

Clone the repository and open the Visual Studio solution file "nrd2dat2.sln".  Build the solution to generate the executable.

### Who do I talk to? ###

richard.gardner@ntnu.no
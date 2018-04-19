# README

nrd2dat is a utility for converting Neuralynx raw data ".nrd" files to the flat binary ".dat" format used by spike sorting packages such as Klusta, Kilosort and Phy. AD channel samples, digital IO port values and timestamps are extracted to separate .dat files. All records are error-checked using CRC32 to ensure validity of the data.

### USAGE

```
nrd2dat <raw_data_file> <channel_map_file> <options>

AVAILABLE OPTIONS:

--input_range <value>
  specify input range of output file in microvolts (default 2000)

--lowcut <value>
  specify low-cut frequency of high-pass filter in Hz (default 0.1)

--buffer_size <value>
 specify size of IO buffers in bytes (default 65536)

--debug
  display debug-level log messsages in console.

--help
  display this text.
```

A Neuralynx acquisition system will save all of its AD channels to the .nrd file, regardless of which channels are actually used. nrd2dat extracts a specified subset of the AD channels, using a channel map indicating the numbers of the required channels. Samples from these channels will be written to the .dat file in the same order they are specified in the channel map.  The channel map must be a text file with the number of one AD channel on each line. Any lines that are empty or  prefixed with '%' will be ignored. See below for an example channel map file.

The Neuralynx raw samples are stored with 24-bit precision, with very large dynamic range (+- 132 mV). When converting the AD signals to the 16-bit output format, signal values are first high-pass filtered to remove near-DC signal offsets, using a zero-phase exponential moving average filter (the same as the 'DCO' high-pass filter used in Cheetah). Next, the filtered signal values are multiplied by a fixed scaling factor determined from the "input_range" parameter, which determines the dynamic range of the output 16-bit signals.

### OUTPUT FILES
Three headerless .dat files are generated, with names based on the input raw data file:

  - **`<basename>_samples.dat`** : a stream of 16-bit signed integers ordered first by AD channel (as specified in the channel map) and second by sample.
  
  - **`<basename>_timestamps.dat`** : a stream of 64-bit unsigned integers storing the timestamp of each record. 
  
  - **`<basename>_ttl.dat`** : values of the digital IO ports for each record, encoded as unsigned 32-bit integers.
  
A log file is also created, with the naming pattern **`<basename>_nrd2dat.log`**.

### EXAMPLE CHANNEL MAP 
Below is an example of a text file that can be used as the channel map argument to nrd2dat. The channel map includes the AD channels from tetrodes 1, 2 and 3 on the Neuralynx HS54 headstage:
```
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
```

### How do I get set up? ###

Clone the repository and open the solution file "nrd2dat2.sln" in Microsoft Visual Studio.  Build the solution to generate the executable.

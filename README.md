# README

nrd2dat is a utility for converting Neuralynx raw data ".nrd" files to the flat binary ".dat" format used by spike sorting packages such as Klusta, Kilosort and Phy.  It uses buffered IO streams for fast reading and writing of large data files.  The speed of your hard drive will likely be the limiting factor for performance.  A 300 GB  .nrd file is typically processed in ~1 hour on my (fairly slow) system.  Error-checking is performed to ensure records have been read correctly from the raw data file.

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

The Neuralynx .nrd files include all available channels, given the system's input board configuration.  This converter requires a channel map indicating the numbers of the AD channels to be extracted.  Samples from these channels will be written to the .dat file in the order specified in the channel map.  The channel map must be a text file with the number of a single AD channel on each line. Any lines that are empty or  prefixed with '%' will be ignored.  See below for an example channel map.

The Neuralynx raw samples are stored with 24-bit precision, with very large dynamic range (+- 132 mV). When converting from 24-bit to 16-bit, signal values are first high-pass filtered to remove near-DC signal offsets, using a zero-phase exponential moving average filter. Next, the filtered signal values are multiplied by a fixed scaling factor determined from the "input_range" parameter which sets the dynamic range of the output 16-bit signals.

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

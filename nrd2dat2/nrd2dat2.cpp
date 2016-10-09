// nrd2dat2.cpp
//
/* File converter from Neuralynx raw data format (.nrd) to flat binary .dat file.

This is a simple and fast way of extracting the sample values of a specified set of AD channels and saving
them in the flat 16-bit binary "dat" format used by Klusta, Kilosort and Phy.  Neuralynx AD converters have
a +/- 132 mV range and work to 24-bit precision.  The 24-bit samples are stored as 32-bit integers in the 
raw data file.  For conversion into 16-bit integers, the sample values are cast to floating precision 
and are high-pass-filtered at 100 Hz to remove all large low-frequency fluctuations from the signals. This
greatly reduces the range of the signal, allowing downcasting to 16-bit integers with only small loss of 
precision.

Each raw data record has a 64-bit unsigned integer timestamp value.  A stream of these values is saved 
separately in a decicated .dat file.

Rich 2016-08-11

This program must be called with two arguments:

1) Name of the source .nrd file

2) Name of a text file containing the channel map.
The channel map must be a text file with the number of a single AD channel on each line. Samples will be 
written to the .dat file in the order specified in the channel map.  Any lines that are empty or  prefixed 
with '%' will be ignored.  See below for an example channel map

The output is a pair of .dat files whose names begin with the base name of the input raw data file.
Both files are headerless.

The first file <basename>_samples.dat contains a stream of 16-bit signed integers ordered first by AD
channel (as specified in the channel map) and second by sample.

The second file <basename>_timestamps.dat contains a stream of 64-bit unsigned integers, representing
the timestamp of each record.


EXAMPLE CHANNEL MAP (corresponding to tetrodes 1, 2 and 3 on the HS54 headstage):

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


#include "stdafx.h"

#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <cmath>
#include <iostream>
#include <time.h>
#include "buffers.h"
#include "filters.h"


#define SAMPLE_SCALE_FACTOR 0.5f  // samples are multipled by this before conversion to int16, to prevent clipping
#define RESONANCE 1.414213562373f  // pi = lowest possible resonance
#define SAMPLE_RATE 32000		   // necessary for calculating filter params
#define HPF_CUTOFF 100.0f		   // set cutoff for HPF applied to samples

#define BUFFER_SIZE 65536		  // set size of buffer used by IO streams
#define RECORD_LIMIT 0			  // limit processing to a fixed number of records
#define RECORD_STX 2048			  // value of STX field at position 0 in each record
#define RECORD_SIZE 584			  // record size in bytes
#define NUM_CHANNELS_TOTAL 128	  // number of channels saved to NRD file
#define HEADER_FIELDS 17		  // number of header (pre-data-packet) fields in each record


/*
FUNCTION DECLARATIONS
*/

// Read signed int32 from specified position in buffer
__int32 getInt32(char* buff, int idx);

// File name construction funcs
string makeOutputBaseFileName(string& inputFileName);
string makeOutputSamplesFileName(string& inputFileName);
string makeOutputTimestampsFileName(string& inputFileName);

// Channel map reader
vector<int> readChannelMap(char* fileName);

bool isNumber(const std::string& s);

// (For debugging)
void printBits(unsigned char byte);
void printBits(unsigned __int64 val);

// Find the next 'STX' value in the file
bool seekNextStx();

// Process one record
bool processRecord(char* sampleBuff, char* writeBuff, int* bufferStartIndices, int& numChannelsToRead, vector<FilterButterworth>& filters);

// Cyclic redundancy error check
__int32 CRC(char* buff);

// Filter one data sample
void filterSample(__int32 val, FilterButterworth& filter, char* buff, int buffIndex);

// Some global vars
int* numClippedSamples;
int currentChannel;
unsigned __int64 lastTimestamp = 0;

char writeBufferTimestamp[8];

// Declare read/write streams as globals
ifstream readStream;
ofstream writeStreamSamples;
ofstream writeStreamTimestamps;

using namespace std;

int main(int argc, char** argv) {

	// Start timer
	clock_t startTime = clock();

	// Check inputs
	if (argc != 3) {
		cout << "Requires raw data and channel map file names." << endl;
		return 0;
	}

	int nChunks = 0;

	cout << endl;

	string inputDataFileName = string(argv[1]);
	string outputFileName = makeOutputSamplesFileName(inputDataFileName);
	string outputTimestampFileName = makeOutputTimestampsFileName(inputDataFileName);

	writeStreamSamples = ofstream(outputFileName.c_str(), ios::binary);
	char bufferSamples[BUFFER_SIZE];
	writeStreamSamples.rdbuf()->pubsetbuf(bufferSamples, BUFFER_SIZE);

	writeStreamTimestamps = ofstream(outputTimestampFileName.c_str(), ios::binary);
	char bufferTimestamps[BUFFER_SIZE];
	writeStreamTimestamps.rdbuf()->pubsetbuf(bufferTimestamps, BUFFER_SIZE);

	// argv[2] should be the channel map file name
	cout << endl;
	vector<int> channelMap = readChannelMap(argv[2]);
	cout << endl;

	// Check that channel map was valid
	if (channelMap.size() == 0) {
		cout << "Invalid channel map. Aborting." << endl;
		return 0;
	}

	// Determine how many channels to read from NRD file
	int numChannels = channelMap.size();

	// Initialize clipped channels array
	numClippedSamples = new int[numChannels];
	for (int ch = 0; ch < numChannels; ch++) {
		numClippedSamples[ch] = 0;
	}

	// argv[1] should be the raw data file name
	cout << "Opening " << argv[1] << endl;
	readStream = ifstream(argv[1], ios::binary);

	// Set the buffer
	char bufferInput[BUFFER_SIZE];
	readStream.rdbuf()->pubsetbuf(bufferInput, BUFFER_SIZE);

	if (readStream.is_open()) {
		cout << "Successfully opened input data file " << inputDataFileName << endl;
	}
	else {
		cout << "Aborting: failed to open raw data file " << argv[1] << endl;
		return 0;
	}

	if (writeStreamSamples.is_open()) {
		cout << "Successfully opened output samples file " << outputFileName << endl;
	}
	else {
		cout << "Aborting: failed to open output samples file " << outputFileName << endl;
		return 0;
	}

	if (writeStreamTimestamps.is_open()) {
		cout << "Successfully opened output timestamps file " << outputTimestampFileName << endl;
	}
	else {
		cout << "Aborting: failed to open output timestamps file " << outputTimestampFileName << endl;
		return 0;
	}

	// Calculate file size
	readStream.seekg(0, readStream.end);
	streamoff fileSize = readStream.tellg();
	readStream.seekg(0, readStream.beg);

	// Create byte buffer with length in terms of uint32
	//cout << "Current input file pos = " << readStream.tellg()  << "." << endl;
	cout << "Scanning for data start" << endl;

	cout << "Current input file pos = " << readStream.tellg() << "." << endl;
	bool foundDataStart = seekNextStx();

	// Now we've reached the start of the data, keep getting records one by one until we reach the end of the file
	char recordBuffer[RECORD_SIZE];

	cout << "Allocating write buffer for " << numChannels << " channels x int16." << endl << endl;
	char* writeBuffer = new char[numChannels*sizeof(__int16)];

	int recordCounter = 0;
	int badRecordCounter = 0;

	bool validRecord;

	// Create an array of start indices for data samples
	cout << "Expected AD sample locations in raw data record:" << endl;
	int* sampleStartIndex = new int[numChannels];
	for (int i = 0; i < numChannels; i++) {
		sampleStartIndex[i] = (HEADER_FIELDS + channelMap.at(i)) * sizeof(__int32);
		printf("AD %d\t @ byte %d\n", channelMap.at(i), sampleStartIndex[i]);
	}

	int percentComplete = 0;
	int checkInterval = 100;

	cout << endl << endl << "Processing raw data records..." << endl;

	// Create filters
	vector<FilterButterworth> filters;
	for (int c = 0; c < numChannels; c++) {
		filters.push_back(FilterButterworth(HPF_CUTOFF, SAMPLE_RATE, RESONANCE));
	}

	bool timeEstimatePrinted = false;

	clock_t recordReadStartTime = clock();

	while (readStream.read(recordBuffer, RECORD_SIZE)) {


		validRecord = processRecord(recordBuffer, writeBuffer, sampleStartIndex, numChannels, filters);

		// Every 100th record, check progress and post to console
		if (recordCounter % checkInterval == 0)  {

			int currentPercentComplete = (float)readStream.tellg() / (float)fileSize * 100.0;

			if (currentPercentComplete > (percentComplete+1)) {

				if (!timeEstimatePrinted) {
					clock_t estimatedTotalTime = 50.0 * (clock() - recordReadStartTime) / CLOCKS_PER_SEC;
					if (estimatedTotalTime < 60) {
						cout << "Estimated total time = " << (float)estimatedTotalTime << " seconds.";
					}
					else {
						cout << "Estimated total time = " << (float)estimatedTotalTime/60 << " minutes.";
					}
					cout << endl << endl;

					// Print progress bar markers
					cout << endl << "0%";
					for (int i = 0; i < 50; i++) {
						cout << "-";
					}
					cout << "100%" << endl << endl;

					cout << "  ";
					timeEstimatePrinted = true;
				}

				// Fill up progress bar
				percentComplete += 2;
				cout << "=";

			}
		}

		if (validRecord) {
			writeStreamSamples.write(writeBuffer, numChannels*sizeof(__int16));
			writeStreamTimestamps.write(writeBufferTimestamp, sizeof(unsigned __int64));
			recordCounter++;
		}
		else {

			badRecordCounter++;
			//debug = true;

			// Seek the next STX value
			//cout << "Seeking next valid record from index " << readBuffer.getIndex() << endl;
			bool foundNext = seekNextStx();
			//cout << "Buffer index after finding record start = " << readBuffer.getIndex() << "." << endl;

			// If the next STX wasn't found, we're done with the file.
			if (!foundNext) {
				break;
			}

		}


		if (RECORD_LIMIT != 0 && recordCounter >= RECORD_LIMIT) {
			break;
		}

	}

	cout << endl << endl;

	// Create new buffer with size matched to a single record
	if (readStream.eof())
	{
		cout << "Finished! " << recordCounter << " samples read, " << badRecordCounter << " bad records.";
		cout << endl << endl;
	}
	else if (readStream.bad())
	{
		// Error reading
		cout << "Error while reading data" << endl << endl;
	}

	// CLOSE FILES
	readStream.close();
	writeStreamSamples.close();
	writeStreamTimestamps.close();

	// Print some summary text
	float meanPercentClipped = 0;
	cout << "Clipped samples:" << endl;
	for (int ch = 0; ch < numChannels; ch++) {
		float percentClipped = (float)numClippedSamples[ch] / (float)recordCounter * 100.0;
		meanPercentClipped += percentClipped/numChannels;
		printf("AD ch %d: \t %d\t (%.3f%%)\n", channelMap[ch], numClippedSamples[ch], percentClipped);
	}

	cout << endl;
	printf("Mean clipping across all AD channels was %.3f%%.\n", meanPercentClipped);
	cout << endl << endl;

	// Get time elapsed
	clock_t endTime = clock();
	clock_t timeElapsed = endTime - startTime;

	cout << "That took " << (float)timeElapsed/CLOCKS_PER_SEC << " seconds." << endl << endl;

	// COLLECT GARBAGE
	delete[] writeBuffer;
	delete[] numClippedSamples;
	delete[] sampleStartIndex;



	return 0;

}

__int32 getInt32(char* buff, int idx) {

	// Reinterpret all buffer bytes as unsigned chars
	unsigned char* bytes = reinterpret_cast<unsigned char*>(buff);

	return (bytes[idx + 3] << 24) + (bytes[idx + 2] << 16) + (bytes[idx + 1] << 8) + (bytes[idx]);

}

unsigned __int32 getUnsignedInt32(char* buff, int idx) {

	// Reinterpret all buffer bytes as unsigned chars
	unsigned char* bytes = reinterpret_cast<unsigned char*>(buff);

	return (bytes[idx + 3] << 24) | (bytes[idx + 2] << 16) | (bytes[idx + 1] << 8) | (bytes[idx]);

}

void printBits(unsigned char byte) {

	char bit;
	cout << "Char value " << (int)byte << ", bits:";

	for (int i = 0; i < 8; i++) {
		bit = (byte >> (7-i)) & 1;
		cout << (int) bit;
	}
	cout << endl;
}

void printBits(unsigned __int64 val) {

	char bit;
	cout << "Char value " << val << ", bits:";

	for (int i = 0; i < 64; i++) {
		bit = (val >> (63-i)) & 1;
		cout << (int)bit;
	}
	cout << endl;
}

bool processRecord(char* sampleBuff, char* writeBuff, int* writeBuffStartIndices, int& numChannelsToRead, vector<FilterButterworth>& filters) {
	// Get sample and timestamp data from a single raw data record

	// Get the packet ID. Reject record if value is not 1
	__int32 packetId = getInt32(sampleBuff, 1 * sizeof(__int32));

	if (packetId != 1) {
		cout << "Invalid packet id" << endl;
		return false;
	}

	// Get the packet size; 
	__int32 packetSize = getInt32(sampleBuff, 2 * sizeof(__int32));

	if (packetSize != NUM_CHANNELS_TOTAL+10) {
		cout << "Invalid packet size of " << packetSize << endl;
		for (int i = 0; i < 4; i++) {
			printBits((unsigned char)sampleBuff[2 * sizeof(__int32) + i]);
		}
		return false;
	}

	// Run CRC
	__int32 crcValue = CRC(sampleBuff);

	if (crcValue != 0) {
		cout << endl << "CRC failed with value of " << crcValue << endl << ".";
		return false;
	}

	// OK! Now we can read the data!

	// Read the timestamp first; it comes in 2 int32 pieces which must
	// be combined into an int64
	unsigned __int32 timestampHigh =  getUnsignedInt32(sampleBuff, 3 * sizeof(__int32));
	unsigned __int32 timestampLow = getUnsignedInt32(sampleBuff, 4 * sizeof(__int32));

	unsigned __int64 timestamp = timestampHigh;
	timestamp <<= 32;
	timestamp += timestampLow;

	int idx;
	for (int i = 0; i < sizeof(unsigned __int64); i++) {
		idx = 3 * sizeof(__int32) + i;
	}


	// Check the timestamp is greater than the previous value.
	// Reject the record if not.
	if (timestamp <= lastTimestamp) {
		cout << "Invalid timestamp of " << timestamp << "." << endl;
		return false;
	}
	else {
		// Fill the timestamp write buffer
		for (int i = 0; i < 8; i++) {
			writeBufferTimestamp[i] = ((timestamp >> (i * 8)));
		}
		
	}

	int numBytes = sizeof(__int32);
	int c = 0;

	__int32 val;

	for (int ch = 0; ch < numChannelsToRead; ch++) {

		// Set global current channel var
		currentChannel = ch;

		// Get the uint32 value of the current sample
		val = getInt32(sampleBuff, writeBuffStartIndices[ch]);

		// Filter sample and enter info write buffer
		filterSample(val, filters.at(ch), writeBuff, c);

		c++;
		c++;

	}

	return true;

}

void filterSample(__int32 val, FilterButterworth& filter, char* buff, int buffIndex) {
	
	filter.Update((float)val);
	
	float newVal = filter.Value() * SAMPLE_SCALE_FACTOR;

	// Clip value if outside of int16 range.
	if (newVal > 32767) {
		newVal = 32767;
		numClippedSamples[currentChannel] += 1;
	}
	else if (newVal < -32768) {
		newVal = -32768;
		numClippedSamples[currentChannel] += 1;
	}
	
	__int16 newValInt = newVal;

	buff[buffIndex] = newValInt;
	buff[buffIndex+1] = newValInt >> 8;
	
}

__int32 CRC(char* buff) {

	int recordHeaderFooterSize = 18;  // in int32
	int recordFieldCount = recordHeaderFooterSize + NUM_CHANNELS_TOTAL;
	int fieldIndex = 0;

	__int32 crcValue = 0;

	__int32 currentField;

	//loop through each field in the record
	for (fieldIndex = 0; fieldIndex<recordFieldCount; fieldIndex++) {

		currentField = getInt32(buff, fieldIndex * sizeof(__int32));
		//logically or all values in packet
		crcValue ^= currentField;
	}

	return crcValue;

}

string makeOutputSamplesFileName(string& inputFileName)
{
	string outputName = makeOutputBaseFileName(inputFileName);
	outputName.append("_samples.dat");
	cout << "Output samples file name = " << (outputName) << endl;
	return outputName;
}

string makeOutputTimestampsFileName(string& inputFileName)
{
	string outputName = makeOutputBaseFileName(inputFileName);
	outputName.append("_timestamps.dat");
	cout << "Output timestamps file name = " << (outputName) << endl;
	return outputName;
}

string makeOutputBaseFileName(string& str) {

	//string str(inputFileName);

	int idx0 = str.rfind("\\");
	int idx1 = str.rfind(".");

	string path = str.substr(0, idx0 + 1);
	string name = str.substr(idx0 + 1, idx1 - idx0 - 1);
	string ext = str.substr(idx1);

	return name;

}

vector<int> readChannelMap(char* fileName) {

	vector<int> channels;

	string line;

	ifstream myfile(fileName);

	int lineCount = 0;

	cout << "Parsing channel map file " << fileName << endl;

	//int currentChannel;

	if (myfile.is_open())
	{
		while (!myfile.eof())
		{
			getline(myfile, line);

			// Skip any blank or comment lines
			if (line.length() == 0 || line.at(0) == '%') {
				continue;
			}

			if (isNumber(line)) {
				channels.push_back(atoi(line.c_str()));
				cout << channels.back() << " ";
				lineCount++;
			}
			else {
				cout << "Invalid line \"" << line << " in channel map file." << endl;
				channels.clear();
				break;
			}


		}
		cout << endl;
		cout << "Read " << lineCount << " lines (channels)." << endl;
		myfile.close();
	}

	else { cout << "Unable to open file"; }

	return channels;

}

bool isNumber(const string& s)
{
	string::const_iterator it = s.begin();
	while (it != s.end() && isdigit(*it)) ++it;
	return !s.empty() && it == s.end();

}

bool seekNextStx() {

	char buffer[sizeof(__int32)];

	int buffIdx = 0;

	while (readStream.read(buffer, sizeof(__int32))) {

		buffIdx++;
		if (getInt32(buffer, 0) == RECORD_STX) {
			cout << "Record start found at byte " << readStream.tellg() << "." << endl;
			// Move back to start of STX
			readStream.seekg(-4, readStream.cur);
			return true;
		}
		// Shouldn't need to seek backwards since all data in file are stored as __int32.
		// We just carry on reading out __int32 until we reach a valid STX value
		//else {
		//	// Move back three, such that next read commences 1 byte after
		//	readStream.seekg(-3, readStream.cur);
		//}
	}

	return false;

}

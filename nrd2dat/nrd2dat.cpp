// nrd2dat.cpp

#include "stdafx.h"

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <cmath>
#include <iostream>
#include <time.h>
#include <algorithm>
#include "buffers.h"
#include "filters.h"
#include "spdlog\spdlog.h"

using namespace std;
namespace spd = spdlog;

// DSP consts
#define RESONANCE 1.414213562373f  // pi = lowest possible resonance
#define SAMPLE_RATE 32000.0        // necessary for calculating filter params
#define HPF_CUTOFF 0.1f            // set low-cut frequency for DCO filter
#define INPUT_RANGE_UV 2000.f      // signal dynamic range
#define FILTER_ORDER 4

// File format consts
#define NBYTES_HEADER 16384        // standard 16 kb Neuralynx file header
#define BUFFER_SIZE 65536		   // default size of buffer used by IO streams
#define RECORD_LIMIT 0			   // limit processing to a fixed number of records
#define RECORD_STX 2048			   // value of STX field at position 0 in each record
#define RECORD_SIZE 584			   // record size in bytes
#define NUM_CHANNELS_TOTAL 128	   // number of channels saved to NRD file
#define HEADER_FIELDS 17		   // number of header (pre-data-packet) fields in each record
#define RAW_AD_BIT_VOLTS 0.132 / 8388608. // +/-132 mV with 24-bit resolution
#define VERSION "0.2.0"
#define VERSION_DATE "06-04-2018"

#ifdef _WIN32
	#define SEP "\\"
#else
    #define SEP "/"
#endif


#define FILTER_TYPE ExpMovingAverage

bool openFiles(string& iName, ifstream& iStream, streamoff& iSize, char*iBuff, ofstream oStreams[], char** oBuffs, size_t bufferSize);

template <typename T> bool getRawHeaderInfo(ifstream& stream, string fieldName, T& val, const T defaultVal);

template <typename T> T string_to_type(const std::string &str);

// Argument parsing
bool parseArgs(int argc, char** argv, bool& debugMode, bool& helpFlag, float& hpfLowCutFreq, float& inputRangeUv, size_t& bufferSize);
template <typename T> T getNumericArg(const char* arg, const string& prmName, T min, T max, const char* units);

// Read signed int32 from specified position in buffer
__int32 getInt32(char* buff, int idx);

// File name construction funcs
string makeOutputFilePath(const string& basePath, const string& fileType, const string& ext);
string makeOutputBaseFilePath(const string& inputFilePath);

// Channel map reader
vector<int> readChannelMap(char* fileName);

// Returns true if a string represents a scalar number
bool isNumber(const std::string& s);

// (For debugging)
void printBits(unsigned char byte);
void printBits(unsigned __int64 val);

// Find the next 'STX' value in the file which marks the beginning of a record
bool seekNextStx(ifstream& readStream);

// Process one record. This is the workhorse function for reading data from the .nrd file.
bool processRecord(char* sampleReadBuff, char* sampleWriteBuff, char* timestampBuff, char* ttlBuff, int* bufferStartIndices, const int& numChannelsToRead, vector<FILTER_TYPE*> filters, int* nClippedSamples, unsigned __int64& lastTimestamp, float& inputRange, size_t& nChannelsTotal);

// Cyclic redundancy error check
__int32 CRC(char* buff, size_t& nChannelsTotal);

float mod_(float x, float div);

// Filter one data sample
void filterSample(float val, FILTER_TYPE* filter, char* buff, int buffIndex, int currentChannel, int* nClippedSamples, float& scaleFactor);


// Logging
template <typename... Args> void lg(const char& level, const char* fmt, const Args &... args);
bool initLogging(const string& inFileName);
std::shared_ptr<spd::logger> loggerConsole;
std::shared_ptr<spd::logger> loggerFile;

string clockTimeToStr(clock_t t);
string basePath = "";
bool loggingInitialized = false;

void printHelp();

void cleanup(ofstream* outFileStreams, ifstream* inFileStream, char* writeBuffer, char* inFileBuffer, int* sampleStartIndex, int* nClippedSamples, char** outFileBuffers);

int main(int argc, char** argv) {

	cout << endl;

	if (argc == 1) {
		cout << "Invalid usage. Refer to documentation below:" << endl << endl;
		printHelp();
		return 0;
	}

	string inFileName;
	try {
		inFileName = string(argv[1]);
		size_t i = inFileName.rfind(SEP, inFileName.length());
		if (i != inFileName.npos) {
			basePath = inFileName.substr(0, i+1);
		}
	} catch (exception& e) {
		cout << e.what();
		cout << "Invalid usage. Refer to documentation below:" << endl << endl;
		printHelp();
		return 0;
	}

	// Check inputs
	if (argc < 3 || argv[1] == "--help") {
		printHelp();
		return 0;
	}

	bool loggingInitialized = initLogging(inFileName);
	if (!loggingInitialized) { return 0; }

	// Parse optional arguments
	float hpfLowCutFreq = HPF_CUTOFF;
	float inputRangeUv = INPUT_RANGE_UV;
	bool debugMode, helpFlag = false;
	size_t bufferSize = BUFFER_SIZE;

	try {
		if (!parseArgs(argc, argv, debugMode, helpFlag, hpfLowCutFreq, inputRangeUv, bufferSize)) { return 0; }
	} catch (exception& e) {
		lg('e', "Parsing args failed: {}. Aborting.", e.what());
		return 0;
	}
	lg('d', "IO buffer size = {} bytes.", bufferSize);

	// Open files
	ifstream* inFileStream = new ifstream();
	ofstream outFileStreams[3];
	streamoff inFileSize;
	char** outFileBuffers = new char*[3]();
	char* inFileBuffer = nullptr;
	char* rawFileHeader = nullptr;
	size_t recordSize, nChannelsTotal;

	if (!openFiles(inFileName, *inFileStream, inFileSize, inFileBuffer, outFileStreams, outFileBuffers, bufferSize)) { 
		cleanup(outFileStreams, inFileStream, nullptr, inFileBuffer, nullptr, nullptr, outFileBuffers);
		return 0;
	};

	string fileVersion;
	double adBitVolts;
	float sampleFrequency;
	getRawHeaderInfo<size_t>(*inFileStream, "-RecordSize", recordSize, RECORD_SIZE);
	getRawHeaderInfo<size_t>(*inFileStream, "-NumADChannels", nChannelsTotal, NUM_CHANNELS_TOTAL);
	getRawHeaderInfo<string>(*inFileStream, "-FileVersion", fileVersion, "(unknown)");
	getRawHeaderInfo<double>(*inFileStream, "-ADBitVolts", adBitVolts, RAW_AD_BIT_VOLTS);
	getRawHeaderInfo<float>(*inFileStream, "-SamplingFrequency", sampleFrequency, SAMPLE_RATE);

	// argv[2] should be the channel map file name
	vector<int> channelMap = readChannelMap(argv[2]);

	// Check that channel map was valid
	const int nChannels = channelMap.size();
	if (nChannels == 0) {
		lg('e', "Invalid channel map. Aborting.");
		cleanup(outFileStreams, inFileStream, nullptr, inFileBuffer, nullptr, nullptr, outFileBuffers);
		return 0;
	}

	// Initialize clipped channels array
	int* nClippedSamples = new int[nChannels];
	for (int ch = 0; ch < nChannels; ch++) { nClippedSamples[ch] = 0; }

	// Scan forward to data start
	lg('d', "Scanning for data start...");
	bool foundDataStart = seekNextStx(*inFileStream);
	lg('d', "Data start found at position {}", inFileStream->tellg());

	// Allocate buffers
	char* recordBuffer = new char[recordSize]; // read buffer for one nrd record
	lg('d', "Allocating write buffer for {} channels x int16.", nChannels);
	char* writeBuffer = new char[nChannels*sizeof(__int16)];
	char writeBufferTimestamp[8];
	char writeBufferTtl[4];

	int recordCounter = 0;
	int badRecordCounter = 0;
	bool validRecord;
	unsigned __int64 lastTimestamp = 0;

	// Create an array of start indices for data samples
	lg('d', "Expected AD sample locations in raw data record:");
	int* sampleStartIndex = new int[nChannels];
	for (int i = 0; i < nChannels; i++) {
		sampleStartIndex[i] = (HEADER_FIELDS + channelMap.at(i)) * sizeof(__int32);
		lg('d', "AD {}\t @ byte {}", channelMap.at(i), sampleStartIndex[i]);
	}

	vector<FILTER_TYPE*> filters;
	for (int c = 0; c < nChannels; c++) {
		FILTER_TYPE* filt = new FILTER_TYPE(hpfLowCutFreq / (sampleFrequency / 2.0));
		filters.push_back(filt);
	}

	double adBitVoltsOut = inputRangeUv / 32767.f / 1e6;
	float scaleFactor = (float)adBitVolts / adBitVoltsOut;
	lg('i', "Input range = {:.0f} uV, scale factor = {:.3f}, bit-volts IN={:.6g}, OUT={:.6g}.", inputRangeUv, scaleFactor, adBitVolts, adBitVoltsOut);
	lg('i', "HPF low-cut frequency = {:.2f} Hz.", hpfLowCutFreq);

	// Prepare to start processing records
	lg('i', "Processing raw data records...");
	clock_t startTime = clock();
	bool timeEstimatePrinted = false;
	int percentComplete = 0;
	int checkInterval = 100;
	bool showProgress = true;

	while (inFileStream->read(recordBuffer, recordSize)) {
		validRecord = processRecord(recordBuffer, writeBuffer, writeBufferTimestamp, writeBufferTtl, sampleStartIndex, nChannels, filters, nClippedSamples, lastTimestamp, scaleFactor, nChannelsTotal);

		// Every 100th record, check progress and post to console
		if (showProgress && recordCounter % checkInterval == 0)  {
			int currentPercentComplete = (float)inFileStream->tellg() / (float)inFileSize * 100.0;
			if (currentPercentComplete > (percentComplete + 1)) {
				if (!timeEstimatePrinted) {
					clock_t estimatedTotalTime = 50.0 * (clock() - startTime);
					lg('i', "Estimated total time {}", clockTimeToStr(estimatedTotalTime));

					// Print progress bar markers
					cout << endl << "0%";
					for (int i = 0; i < 50; i++) { cout << "-"; }
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
			outFileStreams[0].write(writeBuffer, nChannels*sizeof(__int16));
			outFileStreams[1].write(writeBufferTimestamp, sizeof(unsigned __int64));
			outFileStreams[2].write(writeBufferTtl, sizeof(unsigned __int32));
			recordCounter++;
		}
		else {
			// Bad record: seek the next STX value.
			// If we don't find another STX, then we're finished with the file.
			lg('e', "Record ID {} was not successfully read. Searching for next valid record...", recordCounter);
			badRecordCounter++;
			showProgress = false;

			// Search for the next STX, indicating the start of another record.
			// If not found, warn user and finish reading.
			if (!seekNextStx(*inFileStream))  {
				lg('w', "Failed to find any more valid records in the file.");
				break;
			}
		}
		if (RECORD_LIMIT != 0 && recordCounter >= RECORD_LIMIT) { break; }
	}

	if (inFileStream->eof()) {
		// Finished successful read
		lg('i', "Finished!  {} samples read, {} bad records.", recordCounter, badRecordCounter);
	} else if (inFileStream->bad()) {
		// Error reading
		lg('e', "EOF not found in raw data file. Check that the file is not corrupted.");
	}


	// Print some summary text
	if (recordCounter > 0) {
		float meanPercentClipped = 0;
		lg('i', "Showing clipped samples for all channels:");
		for (int ch = 0; ch < nChannels; ch++) {
			float percentClipped = (float)nClippedSamples[ch] / (float)recordCounter * 100.0;
			meanPercentClipped += percentClipped / nChannels;
			lg('i', "AD ch {}: \t {}\t ({:.3f}%)", channelMap[ch], nClippedSamples[ch], percentClipped);
		}
		lg('i', "Mean clipping across all AD channels was {:.3f}%.", meanPercentClipped);
	}

	// Get time elapsed
	clock_t timeElapsed = clock() - startTime;
	lg('i', "That took {}.", clockTimeToStr(timeElapsed));

	// Clean up
	cleanup(outFileStreams, inFileStream, writeBuffer, inFileBuffer, sampleStartIndex, nClippedSamples, outFileBuffers);

	return 0;
}

void cleanup(ofstream* outFileStreams, ifstream* inFileStream, char* writeBuffer, char* inFileBuffer, int* sampleStartIndex, int* nClippedSamples, char** outFileBuffers) {
	// Clean up
	if (inFileStream != nullptr) {
		inFileStream->close();
	}
	if (outFileStreams != nullptr) {
		for (int i = 0; i < 3; i++) { outFileStreams[i].close(); }
	}

	delete[] writeBuffer, inFileBuffer, nClippedSamples, sampleStartIndex, outFileStreams;

	if (outFileBuffers != nullptr) {
		for (int i = 0; i < 3; i++) {
			if (outFileBuffers[i] != nullptr) {
				delete outFileBuffers[i];
			}
		}
	}
	spd::drop_all();
}

string clockTimeToStr(clock_t t) {
	float tSec = (float) t / CLOCKS_PER_SEC;
	char buff[100];
	if (tSec < 60) {
		sprintf_s(buff, "%u seconds", (int) floor(tSec));
	}
	else if (tSec < 3600) {
		sprintf_s(buff, "%u minutes, %u seconds", (int) floor(tSec / 60), (int) mod_(tSec, 60));
	} else {
		sprintf_s(buff, "%u hours, %u minutes, %u seconds", (int) floor(tSec / 3600), (int) floor(mod_(tSec, 3600)/60), (int) floor(mod_(tSec, 60)));
	}
	return buff;
}

float mod_(float x, float div) {
	// Same form used in MATLAB
	float i = 0;
	return div*modf(x / div, &i);
}


__int32 getInt32(char* buff, int idx) {
	// Reinterpret all buffer bytes as signed chars
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
		bit = (byte >> (7 - i)) & 1;
		cout << (int)bit;
	}
	cout << endl;
}


void printBits(unsigned __int64 val) {
	char bit;
	cout << "Char value " << val << ", bits:";
	for (int i = 0; i < 64; i++) {
		bit = (val >> (63 - i)) & 1;
		cout << (int)bit;
	}
	cout << endl;
}


bool processRecord(char* sampleReadBuff, char* writeBuff, char* timestampBuff, char* ttlBuff, int* writeBuffStartIndices, const int& numChannelsToRead, vector<FILTER_TYPE*> filters, int* nClippedSamples, unsigned __int64& lastTimestamp, float& scaleFactor, size_t& nChannelsTotal) {
	// Get sample and timestamp data from a single raw data record

	// Get the packet ID. Reject record if value is not 1
	__int32 packetId = getInt32(sampleReadBuff, 1 * sizeof(__int32));

	if (packetId != 1) {
		lg('e', "Invalid packet id vaue \"{}\".", packetId);
		return false;
	}

	// Get the packet size; 
	__int32 packetSize = getInt32(sampleReadBuff, 2 * sizeof(__int32));

	if (packetSize != nChannelsTotal + 10) {
		lg('e', "Invalid packet size of {}.", packetSize);
		for (int i = 0; i < 4; i++) {
			printBits((unsigned char)sampleReadBuff[2 * sizeof(__int32) + i]);
		}
		return false;
	}

	// Run CRC
	__int32 crcValue = CRC(sampleReadBuff, nChannelsTotal);

	if (crcValue != 0) {
		cout << endl;
		lg('e', "CRC failed with value of {}. Continuing to read remaining records...", crcValue);
		return false;
	}

	// Read the timestamp first; it comes in 2 int32 pieces which must
	// be combined into an int64
	unsigned __int32 timestampHigh = getUnsignedInt32(sampleReadBuff, 3 * sizeof(__int32));
	unsigned __int32 timestampLow = getUnsignedInt32(sampleReadBuff, 4 * sizeof(__int32));
	unsigned __int64 timestamp = timestampHigh;
	timestamp <<= 32;
	timestamp += timestampLow;

	// Check the timestamp is greater than the previous value.
	// Reject the record if not.
	if (timestamp <= lastTimestamp) {
		lg('e', "Invalid timestamp of {}.", timestamp);
		return false;
	}
	else {
		// Fill the timestamp write buffer
		for (int i = 0; i < 8; i++) {
			timestampBuff[i] = ((timestamp >> (i * 8)));
		}
		lastTimestamp = timestamp;
	}

	// Get the parallel TTL input port value
	unsigned __int32 ttlInput = getUnsignedInt32(sampleReadBuff, 6 * sizeof(__int32));
	for (int i = 0; i < 4; i++) {
		ttlBuff[i] = ((ttlInput >> (i * 8)));
	}

	int numBytes = sizeof(__int32);
	int c = 0;
	__int32 val;

	for (int ch = 0; ch < numChannelsToRead; ch++) {
		// Get the uint32 value of the current sample
		val = (float)getInt32(sampleReadBuff, writeBuffStartIndices[ch]);
		filterSample(val, filters.at(ch), writeBuff, c, ch, nClippedSamples, scaleFactor);
		c += 2;
	}
	return true;
}


void filterSample(float val, FILTER_TYPE* filter, char* buff, int buffIndex, int currentChannel, int* nClippedSamples, float& scaleFactor) {

	filter->update(val);
	float newVal = (val - filter->value()) * scaleFactor;

	// Clip value if outside of int16 range.
	if (newVal > 32767) {
		newVal = 32767;
		nClippedSamples[currentChannel] += 1;
	}
	else if (newVal < -32768) {
		newVal = -32768;
		nClippedSamples[currentChannel] += 1;
	}

	__int16 newValInt = newVal;
	buff[buffIndex] = newValInt;
	buff[buffIndex + 1] = newValInt >> 8;
}


__int32 CRC(char* buff, size_t& nChannelsTotal) {
	int recordHeaderFooterSize = 18;  // in int32
	int recordFieldCount = recordHeaderFooterSize + nChannelsTotal;
	int fieldIndex = 0;
	__int32 crcValue = 0;
	__int32 currentField;
	//loop through each field in the record and apply bitwise XOR
	for (fieldIndex = 0; fieldIndex < recordFieldCount; fieldIndex++) {
		currentField = getInt32(buff, fieldIndex * sizeof(__int32));
		crcValue ^= currentField;
	}
	return crcValue;
}


string makeOutputFilePath(const string& basePath, const string& fileType, const string& ext) {
	string outputPath = basePath;
	outputPath.append("_");
	outputPath.append(fileType);
	outputPath.append(ext);
	lg('d', "Output {} file path = \"{}\".", fileType, outputPath);
	return outputPath;
}


string makeOutputBaseFilePath(const string& inputFilepath) {

	if (inputFilepath.length() == 0) {
		throw exception("Empty source.");
	}

	size_t idx0 = inputFilepath.rfind("\\");
	size_t idx1 = inputFilepath.rfind(".");

	bool hasPath = idx0 != string::npos;
	bool hasExt = idx1 != string::npos;
	string basepath = "";
	lg('d', "Input file has path = {}.", hasPath);

	if (!hasExt) { idx1 = inputFilepath.length(); }
	lg('d', "Input file has ext = {}.", hasExt);

	basepath = inputFilepath.substr(0, idx1);
	lg('d', "Base filename identified as \"{}\".", basepath);
	return basepath;
}


vector<int> readChannelMap(char* fileName) {

	vector<int> channels;
	string line;
	string channelList = "";
	ifstream file(fileName);

	if (!file.is_open()) {
		lg('e', "Unable to open channel map file.");
		return channels;
	}

	int nChannels = 0;
	lg('d', "Parsing channel map file \"{}\".", fileName);

	while (!file.eof()) {
		getline(file, line);
		// Skip any blank or commented lines
		if (line.length() == 0 || line.at(0) == '%') { continue; }

		if (isNumber(line)) {
			channels.push_back(atoi(line.c_str()));
			channelList.append(line);
			channelList.append(", ");
			nChannels++;
		} else {
			lg('e', "Invalid line in channel map file: \"{}\".", line);
			channels.clear();
			break;
		}
	}

	if (!channelList.empty()){
		channelList.erase(channelList.length() - 2, 2);
		channelList.append(".");
		lg('i', "Read {} channels: {}.", nChannels, channelList);
	} else {
		lg('w', "No channels were read.");
	}
	file.close();

	return channels;
}


bool isNumber(const string& s) {
	bool validChars = !s.empty() && s.find_first_not_of("+-0123456789.") == std::string::npos;
	size_t nPt = std::count(s.begin(), s.end(), '.');
	const char signs[2] = {'+', '-'};
	bool validSign = true;
	for (int i = 0; i < 2; i++) {
		size_t idx = s.find_last_of(signs[i]);
		validSign &= (idx == 0 || idx == string::npos);
	}
	return validChars && nPt <= 1;
}


bool seekNextStx(ifstream& readStream) {
	char buffer[sizeof(__int32)];
	int buffIdx = 0;
	while (readStream.read(buffer, sizeof(__int32))) {
		buffIdx++;
		if (getInt32(buffer, 0) == RECORD_STX) {
			lg('d', "Record start found at byte {}.", readStream.tellg());
			// Move back to start of STX
			readStream.seekg(-4, readStream.cur);
			return true;
		}
	}
	return false;
}


template <typename... Args> void lg(const char& level, const char* fmt, const Args &... args) {

	// Don't log if uninitialized
	if (!loggingInitialized){ return; }

	spd::level::level_enum level_;
	switch (level) {
	case 'v':
		level_ = spd::level::trace;
		break;
	case 'd':
		level_ = spd::level::debug;
		break;
	case 'i':
		level_ = spd::level::info;
		break;
	case 'w':
		level_ = spd::level::warn;
		break;
	case 'e':
		level_ = spd::level::err;
		break;
	case 'c':
		level_ = spd::level::critical;
		break;
	}
	loggerConsole->log(level_, fmt, args...);
	loggerFile->log(level_, fmt, args...);
}


void printHelp() {
	cout << endl;
	cout << "===========================================================" << endl;
	cout << "USAGE:" << endl << endl;
	cout << "nrd2dat <raw_data_file> <channel_map_file> <options>" << endl << endl;
	cout << "AVAILABLE OPTIONS:" << endl << endl;
	cout << "--input_range <value>" << endl << "  specify input range of output file in microvolts (default 2000)" << endl << endl;
	cout << "--lowcut <value>" << endl << "  specify low-cut frequency of high-pass filter in Hz (default 0.1)" << endl << endl;
	cout << "--buffer_size <value>" << endl << " specify size of IO buffers in bytes (default 65536)" << endl << endl;
	cout << "--debug" << endl << "  display debug-level log messsages in console." << endl << endl;
	cout << "--help" << endl << "  display this text." << endl;
	cout << "===========================================================" << endl << endl;
}


bool initLogging(const string& inFileName) {
	string basePath = makeOutputBaseFilePath(inFileName);
	const string logFilePath = makeOutputFilePath(basePath, "nrd2dat", ".log");
	try {
		loggerFile = spd::basic_logger_st("logfile", logFilePath, true);
	}
	catch (const spd::spdlog_ex& e){
		cout << "Log init failed: " << e.what() << endl;
		return false;
	}
	loggerFile->set_level(spd::level::debug);
	loggerConsole = spd::stdout_logger_st("console");
	loggerConsole->set_level(spd::level::info);
	spd::set_pattern("[%H:%M:%S] [%L] %v");
	loggingInitialized = true;
	lg('d', "Base path = {}.", basePath);
	lg('i', "nrd2dat v{}, date {}.", VERSION, VERSION_DATE);
	lg('i', "Saving log file in {}.", logFilePath);
	return true;
}


bool openFiles(string& iName, ifstream& iStream, streamoff& iSize, char* iBuff, ofstream oStreams[3], char* oBuffs[3], size_t bufferSize) {

	// Open the input file
	lg('d', "Opening input raw data file \"{}\".", iName);
	iStream.open(iName, ios::binary);

	// Set the read stream buffer
	iBuff = new char[bufferSize];
	iStream.rdbuf()->pubsetbuf(iBuff, bufferSize);
	if (iStream.is_open()) {
		lg('d', "Raw data file successfully opened.");
	} else {
		lg('e', "Failed to open raw data file. Aborting.");
		return false;
	}

	// Calculate file size
	iStream.seekg(0, iStream.end);
	iSize = iStream.tellg();
	iStream.seekg(0, iStream.beg);

	// Generate output filenames
	string baseFilePath;
	string outputFileTypes[3] = { "samples", "timestamps", "ttl" };
	string outputFilePaths[3];

	try {
		baseFilePath = makeOutputBaseFilePath(iName);
	} catch (std::exception& e) {
		lg('e', "Failed to build base filename: \"{}\".", e.what());
		string msg = "Cannot build base filename: ";
		msg.append(e.what());
		throw exception(msg.c_str());;
	}

	for (int i = 0; i < 3; i++) {
		string fileType = outputFileTypes[i];
		string filename;
		try {
			filename = makeOutputFilePath(baseFilePath, fileType, ".dat");
			outputFilePaths[i] = filename;
		} catch (exception& e) {
			lg('e', "Failed to build output {} filename: \"{}\". Aborting.", fileType, e.what());
			return false;
		}

		try {
			oStreams[i].open(filename.c_str(), ios::binary);
		} catch (exception& e) {
			lg('e', "Failed to create output stream for {} file \"{}\". Aborting.", fileType, e.what());
			return false;
		}

		if (oStreams[i].is_open()) {
			lg('d', "Successfully opened output {} file.", fileType);
		} else {
			lg('e', "Aborting: failed to open output {} file.", fileType);
			return false;
		}

		char* buffer = new char[bufferSize];
		oBuffs[i] = buffer;
		oStreams[i].rdbuf()->pubsetbuf(buffer, bufferSize);
	}

	for (int i = 0; i < 3; i++) {
		if (!oStreams[i].is_open()) {
			lg('e', "Output file {} not open!", outputFilePaths[i]);
			return false;
		}
	}
	return true;
}

template <typename T> bool getRawHeaderInfo(ifstream& stream, string fieldName, T& val, const T defaultVal) {
	
	lg('d', "Scanning Nlx raw data file header for field \"{}\"...", fieldName);

	size_t fieldLen = fieldName.length();

	stream.seekg(0, ios::beg);
	char buff[NBYTES_HEADER];
	stream.read(buff, NBYTES_HEADER);
	string header = string(buff);
	stream.seekg(0, ios::beg);

	size_t idxField = header.find(fieldName);
	if (idxField == header.npos) {
		lg('d', "Header field \"{}\" was not found. Using default value of {}.", fieldName, defaultVal);
		val = defaultVal;
		return false;
	}
	else {
		size_t idxNewline = header.find('\r\n', idxField);
		const string headerVal = header.substr(idxField + fieldLen, idxNewline - idxField - fieldLen);
		val = string_to_type<T>(headerVal);
		lg('d', "Header field \"{}\" value read as {}.", fieldName, val);
		return true;
	}
}

template <typename T> T string_to_type(const std::string &str)
{
	std::istringstream ss(str);
	T num;
	ss >> num;
	return num;
}

bool parseArgs(int argc, char** argv, bool& debugMode, bool& helpFlag, float& hpfLowCutFreq, float& inputRangeUv, size_t& bufferSize) {
	int i = 3;
	while (i < argc) {
		string prm = argv[i];

		// No-argument switches
		if (prm == "--debug") {
			loggerConsole->set_level(spd::level::debug);
			i++;
			continue;
		} else if (prm == "--help") {
			printHelp();
			return 0;
		}

		// Other switches require arguments
		if (argc < i + 2) {
			lg('e', "Option \"{}\" must be followed by a corresponding value.", prm);
			return 0;
		}

		const char* val = argv[i + 1];

		if      (prm == "--lowcut")      { hpfLowCutFreq = getNumericArg<float>(val, prm, 0.1, 15000, "Hz"); }
		else if (prm == "--input_range") { inputRangeUv = getNumericArg<float>(val, prm, 1, numeric_limits<float>::max(), "uV"); }
		else if (prm == "--buffer_size") { bufferSize = getNumericArg<int>(val, prm, 1, numeric_limits<int>::max(), "bytes"); }
		else { lg('e', "Invalid optional argument \"{}\".", prm); return 0; }
		i += 2;
	}
	return true;
}

template <typename T> T getNumericArg(const char* arg, const string& prmName, T min, T max, const char* units) {
	T val;
	if (isNumber(arg)) {
		lg('d', "Parameter \"{}\" value \"{}\" numeric=YES.", prmName, arg);
		if (is_integral<T>::value){
			val = atoi(arg);
		} else if(is_floating_point<T>::value){
			val = (float)atof(arg);
		}

		if (val < min || val > max) {
			lg('e', "Parameter \"{}\" value \"{}\" is out of permitted bounds ({} ... {} {})", prmName, arg, min, max, units);
			throw exception("Out-of-bounds argument.");
		}
	}
	else {
		lg('e', "Failed to interpret integer value \"{}\" of parameter \"{}\".", arg, prmName);
		throw exception("Failed to get integer arg.");
	}
	return val;
}
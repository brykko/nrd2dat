#ifndef _BUFFERS_H_INCLUDED_
#define _BUFFERS_H_INCLUDED_

#define BUFFER_LENGTH 65536 // 64 kB
#include <iostream>
#include <fstream>

using namespace std;

class Buffer {


protected:

	char bytes[BUFFER_LENGTH];
	long int index = 0;
	long int length = BUFFER_LENGTH;

	streamoff filePos;

public:

	Buffer(): bytes() {}


	streamoff getPositionInFile() {
		return filePos;
	}

	streamoff getIndex() {
		return index;
	}


	streamoff getPosition() {
		// File position corresponds to the END of the current cycle.
		// Net position must subtract the unread portion of the buffer
		return getPositionInFile() - length + index + 1;
	}

	// Public read method for filling a given buffer
	bool transfer(char buffer[], int nBytes) {

		// Determine if we will need to refill the buffer
		if (nBytes + index > length) {

			//cout << "Refilling!" << endl;

			int nBytes1 = length - index;
			int nBytes2 = nBytes - nBytes1;

			// Get all bytes up to the end of the buffer
			transferNBytes(buffer, nBytes1, 0);

			// Now we've reached the end, refill
			cycle();
			index = 0;

			// Read the rest of the bytes
			transferNBytes(buffer, nBytes2, nBytes1);

		}

		// If no refill needed, just get all
		else {
			//cout << " no refill needed." << endl;
			transferNBytes(buffer, nBytes, 0);
		}

		return true;

	}

protected:

	// Internal cycle method for cycling (i.e. refreshing/purging) the buffer
	virtual void cycle() = 0;

	// Internal method for transferrring a set number of bytes between the internal buffer and a smaller user buffer
	virtual void transferNBytes(char buffer[], int numBytes, int startIndex) = 0;


};


class ReadBuffer : public Buffer {

private:


	streamoff fileLength;
	ifstream* stream;

public:

	ReadBuffer() : Buffer() {
		stream = 0;
	}

	explicit ReadBuffer(ifstream& str) {

		stream = &str;

		cout << "Creating ReadBuffer!" << endl;

		// Get current position and length of file
		filePos = stream->tellg();
		stream->seekg(0, stream->end);
		fileLength = stream->tellg();
		stream->seekg(filePos, stream->beg);

		// Fill the buffer
		cycle();

	}

	// Public read method for filling a given buffer
	bool transfer(char buffer[], int nBytes) {
		Buffer::transfer(buffer, nBytes);
		return getRemaining() > 0;
	}

	// Shift index by a certain number of bytes
	bool shift(int n) {

		// Positive shift
		if (n > 0) {

			// If shift spills over end of file, return false
			if (n + index >= fileLength) {
				return false;
			}

			// Check for spilling over boundary
			if (n >= (length - index)) {
				index -= length;
				cycle();
			}

			// Shift is OK - apply it
			index += n;

		}

		// Negative shift
		else if (n < 0) {

			// If shift goes before file beginning, return false
			if (n + index < 0) {
				return false;
			}

			// Check for spilling over boundary
			if (n + index < 0) {
				index += length;
				reverseCycle();
			}

			// Shift is OK - apply it
			index += n;

		}

		index = index % length;


	}

	streamoff getFileLength() {
		return fileLength;
	}

	// Remaining bytes in file
	streamoff getRemainingInFile() {
		return fileLength - filePos - 1;
	}

	// Total bytes remaining
	streamoff getRemaining() {
		streamoff bytesInFile = getRemainingInFile();
		if (bytesInFile > 0) {
			return bytesInFile + (length - index);
		}
		else {
			return (fileLength % length) - index - 1;
		}
		return getRemainingInFile() + (length - index);
	}


	bool good() {
		return stream->good();
	}

	bool isDone() {
		return (getRemaining() == 0);
	}

private:

	// Get a specified number of bytes from the internal buffer
	void transferNBytes(char buffer[], int numBytes, int startIndex) {
		for (int n = 0; n < numBytes; n++) {
			buffer[startIndex + n] = bytes[index++];
		}
	}

	// Refill the internal buffer
	void cycle() {

		//cout << "Stream OK = " << stream->good() << endl;
		// If number of remaining bytes exceeds or equals length of buffer, fill whole buffer
		if (getRemaining() > length) {
			stream->read(bytes, length);
			filePos += length;
		}
		// Otherwise, get all remaining bytes
		else {
			stream->read(bytes, getRemaining());
			filePos = fileLength - 1;
		}

	}

	void reverseCycle() {
		// Check that the file position is advanced enough to reverse one cycle
		if (filePos >= (length - 1)) {

			// Move back TWO cycles
			stream->seekg(-2 * length, stream->cur);
			filePos -= 2 * length;

			// Now refill from current position
			cycle();

		}
	}


};

class WriteBuffer : public Buffer {

private:

	ofstream* stream;

public:

	WriteBuffer() : Buffer() {
		stream = 0;
	}

	explicit WriteBuffer(ofstream& str) {
		stream = &str;
		// Get current position and length of file
		filePos = stream->tellp();
	}

	// Public transfer method for filling from a given buffer
	bool transfer(char buffer[], int nBytes) {
		Buffer::transfer(buffer, nBytes);
		return true;
	}

	bool good() {
		return stream->good();
	}

private:

	// Put a specified number of bytes into the internal buffer
	void transferNBytes(char buffer[], int numBytes, int startIndex) {
		for (int n = 0; n < numBytes; n++) {
			bytes[index++] = buffer[startIndex + n];
		}
	}

	// Refill the internal buffer
	void cycle() {
		// Write all bytes up to the current index
		stream->write(bytes, index + 1);
		filePos += index;
	}


};

#endif
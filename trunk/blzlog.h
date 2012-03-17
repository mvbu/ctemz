#ifndef _INCL_BLZLOG_H_
#define _INCL_BLZLOG_H_

#include <string>
#include <sstream>
#include <cstring>

using namespace std;

typedef enum {
	PSYCHO,
	DEBUG,
	INFO,
	WARN,
	ERROR
} BlzLogLevel;

static const int LEVEL_COUNT = 5;

/**
	 Simple encapsulation of some logging functionality. Has logging levels and outputs to cout.
 */
class BlzLog {
 private:
	static BlzLogLevel level;
 public:
	// Set the lowest level of logs to display. DEBUG lowest, ERROR highest
	static void setLevel(BlzLogLevel level_);
	static void setLevel(char* level_);
	// So that anyone can make use of the log level in an 'if' statement
	static BlzLogLevel getLevel();
	static std::string getLevelString();
	static void output(BlzLogLevel level_, string msg);
	static void output(BlzLogLevel level_, ostringstream& msg);

	// Convenience methods
	static void debug(string msg);
	static void info(string msg);
	static void warn(string msg);
	static void error(string msg);
	static void debug(ostringstream& msg);

	template<class T>
  static void logArray(string name, T data[], int size, int rowLength, BlzLogLevel level_) {
		int i;
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ":" << endl;
		int iCol = 0;
		for(i=0; i<size;i++) {
			oss << data[i] << " ";
			iCol++;
			if(iCol==rowLength) {
				oss << endl;
				iCol=0;
			}
		}    
		output(level_, oss);
	}

	template<class T>
	static void debugArray(string name, T data[], int size, int rowLength) {
		logArray(DEBUG, name, data, size, rowLength);
	}

	template<class T>
	static void psychoArray(string name, T data[], int size, int rowLength) {
		logArray(PSYCHO, name, data, size, rowLength);
	}

	template<class T>
		static void outputVector(BlzLogLevel level_, string name, T data[], int size) {
		int i;
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ": ";
		for(i=0; i<size;i++) {
			T value = data[i];
			if(abs(value) < 1e-150) value = (T)0;
			oss << value << " ";
		}    
		output(level, oss);
	}

	template<class T>
	static void debugVector(string name, T data[], int size) {
		outputVector(DEBUG, name, data, size);
	}

	template<class T>
	static void outputComplexVector(BlzLogLevel level_, string name, T data[], int size) {
		int i;
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ":" << endl;
		for(i=0; i<size;i++) {
			T realvalue = data[i*2];
			if(abs(realvalue) < 1e-150) realvalue = (T)0;
			T imagvalue = data[i*2+1];
			if(abs(imagvalue) < 1e-150) imagvalue = (T)0;
			oss << realvalue << " " << imagvalue << endl;
		}    
		output(level_, oss);
	}

	template<class T>
	static void debugComplexVector(string name, T data[], int size) {
		outputComplexVector(DEBUG, name, data, size);
	}

	
	template <class T>
  static void outputScalar(BlzLogLevel level_, string name, T scalar) {
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ": ";
		oss << scalar;
		output(level_, oss);
	}	

	template <class T>
  static void psychoScalar(string name, T scalar) {
		outputScalar(PSYCHO, name, scalar);
	}
	template <class T>
  static void debugScalar(string name, T scalar) {
		outputScalar(DEBUG, name, scalar);
	}

	template <class T>
  static void infoScalar(string name, T scalar) {
		outputScalar(INFO, name, scalar);
	}

	template <class T>
  static void errorScalar(string name, T scalar) {
		outputScalar(ERROR, name, scalar);
	}

	template <class T>
	static void debugComplex(string name, T real, T imag) {
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ": ";
		oss << real;
		if(imag != 0.0) {
			if(imag > 0)
				oss << "+";
			oss << imag << "i";
		}
		output(DEBUG, oss);
	}	

	template <class T>
	static void outputScalarPair(BlzLogLevel level_, string name, T scalar1, T scalar2) {
		ostringstream oss;
		oss.precision(10);
		if(!name.empty())
			oss << name << ": ";
		oss << scalar1;
		oss << " ";
		oss << scalar2;
		output(level_, oss);
	}	

	template <class T>
	static void debugScalarPair(string name, T scalar1, T scalar2) {
		outputScalarPair(DEBUG, name, scalar1, scalar2);
	}	

 private:
  BlzLog();
  ~BlzLog();
};
    
#endif // _INCL_BLZLOG_H_


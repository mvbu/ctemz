#include <iostream>
#include <sstream>
#include "blzlog.h"

using namespace std;

static char LEVEL_STRINGS[][32] = {
	"psycho",
	"debug",
	"info",
	"warn",
	"error"
};


BlzLogLevel BlzLog::level = WARN;

BlzLog::BlzLog() {}
BlzLog::~BlzLog() {}

BlzLogLevel BlzLog::getLevel() {
	return level;
}

std::string BlzLog::getLevelString() {
	return std::string(LEVEL_STRINGS[level]);
}

void BlzLog::setLevel(BlzLogLevel level_) {
	level = level_;
}


void BlzLog::setLevel(char* ll) {
	BlzLogLevel level;
	int i = 0;
	for(i=0;i < LEVEL_COUNT;i++) {
		if(strcmp(ll, LEVEL_STRINGS[i])==0) {
			level = (BlzLogLevel)i;
			break;
		}
	}
	
	if(i == LEVEL_COUNT) {
		cout << "% Could not find log level corresponding to " << ll << ". Setting to level 'warning'" << endl;
		// default to WARN level
		level = WARN;
	}

	setLevel(level);
	cout << "% Set logging level to: " << level << " (" << LEVEL_STRINGS[(int)level] << ")" << endl;
}

void BlzLog::output(BlzLogLevel level_, string msg) {
	if(level <= level_)
		cout << msg << endl;
}

void BlzLog::output(BlzLogLevel level_, ostringstream& msg) {
	if(level <= level_)
		cout << msg.str() << endl;
}

void BlzLog::debug(string msg) {
	output(DEBUG, msg);
}

void BlzLog::info(string msg) {
	output(INFO, msg);
}

void BlzLog::warn(string msg) {
	output(WARN, msg);
}

void BlzLog::error(string msg) {
	output(ERROR, msg);
}

void BlzLog::debug(ostringstream& msg) {
	output(DEBUG, msg);
}


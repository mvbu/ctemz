#ifndef _INCL_BLZSIMINPUTREADER_H_
#define _INCL_BLZSIMINPUTREADER_H_

#include <string>
#include "blzsiminput.h"

class BlzSimInputReader {
 public:
  // For now, just one constructor, to set up reading the values from a text file
  BlzSimInputReader(const std::string _filepath);
  // Fill in simInput with the values from (for now) the specified text file
  void read(BlzSimInput& simInput);

 private:
  std::string filepath;
};

#endif // _INCL_BLZSIMINPUTREADER_H_

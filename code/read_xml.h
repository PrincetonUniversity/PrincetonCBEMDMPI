/**
 Read HOOMD XML file format into System object
 \author Nathan A. Mahynski
 **/

#ifndef READ_XML_H_
#define READ_XML_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "system.h"
#include "misc.h"

using namespace std;
using namespace system;
using namespace misc;

/*!
 Parse an XML file to obtain atom information. Returns 0 if successful, -1 if failure.
 \param [in] filename Name of file to open and read.
 */
int System::read_xml (const string filename) {
	FILE *fp1 = mfopen(filename.c_str(), "r");
	;
}

#endif
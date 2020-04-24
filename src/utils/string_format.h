#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Convert double to string with specified number of places after the decimal
   and left padding. */
std::string prd(const double x, const int decDigits, const int width);

std::string center(const std::string s, const int w);
/* Right-aligns string within a field of width w. Pads with blank spaces
   to enforce alignment. */
std::string right(const std::string s, const int w);

/*! Left-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
std::string left(const std::string s, const int w);

#pragma once
#include <stdlib.h>

//////////////////// Utility Functions ////////////////////
inline float RandNext() {
	return static_cast<float>(rand()) / (RAND_MAX + 1);
}
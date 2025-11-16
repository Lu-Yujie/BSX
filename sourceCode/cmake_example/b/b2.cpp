#include "b2.h"
#include "a/a1/a1.h"
#include <iostream>

void b2_function() {
    a1_function();
    std::cout << "This is a1_function() from b2.cpp" << std::endl;
    std::cout << "This is b2_function() from b2.cpp" << std::endl;
}

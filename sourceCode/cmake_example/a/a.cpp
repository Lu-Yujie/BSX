#include "a.h"
#include "a1/a1.h"
#include "a2/a2.h"
#include "a3/a3.h"
#include <iostream>

void a_function() {
    std::cout << "This is a_function() from a.cpp" << std::endl;
    a1_function();
    std::cout << "This is a1_function() from a.cpp" << std::endl;
    a2_function();
    std::cout << "This is a2_function() from a.cpp" << std::endl;
    a3_function();
    std::cout << "This is a3_function() from a.cpp" << std::endl;
}

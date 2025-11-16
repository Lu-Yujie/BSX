#include <iostream>
#include "a/a1/a1.h"
#include "a/a2/a2.h"
#include "a/a3/a3.h"
// #include "a/a.h"
#include "b/b1.h"
#include "b/b2.h"

int main() {
    a1_function();
    a2_function();
    a3_function();
    // a_function();

    b1_function();
    b2_function();

    std::cout << "Main program executed successfully!" << std::endl;
    return 0;
}

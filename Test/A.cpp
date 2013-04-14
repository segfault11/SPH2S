#include <iostream>
struct A 
{
    float Data[2];
};

int main(int argc, const char *argv[])
{
    A a;
    a.Data[0] = 1.0f;
    a.Data[1] = 42.0f;

    A b(a);

    std::cout << b.Data[0] << " " << b.Data[1] << std::endl;


    return 0;
}

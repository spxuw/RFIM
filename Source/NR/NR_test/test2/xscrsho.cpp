#include "nr.h"

// Driver for routine scrsho

DP fx(const DP x)
{
    // return NR::bessj0(x);
    return x*x*x*x - 6*x*x*x + 12.25 *x*x -11*x +3.75;
}

int main(void)
{
        NR::scrsho(fx);
        return 0;
}

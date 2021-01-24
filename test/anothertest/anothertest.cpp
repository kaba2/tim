#include "tim/core/differential_entropy.h"
#include "tim/core/signal_generate.h"

#include <chrono>

using namespace Tim;

template<typename F, typename ...Args>
static auto duration(F&& func, Args&&... args)
{
    using TimeT = std::chrono::seconds;
    auto start = std::chrono::steady_clock::now();
    std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
    return std::chrono::duration_cast<TimeT>(std::chrono::steady_clock::now()-start);
} 

int main() {
    auto signal = generateGaussian(10, 100000);
    SignalData signals[] = { signal };

    auto f = [&]() {
        return differentialEntropyKl(signals);
    };
    
    std::cout << duration(f).count() << std::endl;
    
    return 0;
}
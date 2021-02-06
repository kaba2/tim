#include "tim/core/differential_entropy.h"
#include "tim/core/tsallis_entropy.h"
#include "tim/core/renyi_entropy.h"
#include "tim/core/mutual_information.h"
#include "tim/core/partial_mutual_information.h"
#include "tim/core/partial_transfer_entropy.h"
#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_generate.h"

#include <tbb/task_scheduler_init.h>

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
    tbb::task_scheduler_init init(1);

    auto signal = generateGaussian(10, 1000);
    Signal signals[] = { (Signal)signal };

    auto f = [&]() {
        auto de = differentialEntropyKl(signals);
        auto te = tsallisEntropyLps(signals);
        auto re = renyiEntropyLps(signals);
        auto mu = mutualInformation(signals, signals, 0, 1);
    };
    
    std::cout << duration(f).count() << std::endl;
    
    return 0;
}
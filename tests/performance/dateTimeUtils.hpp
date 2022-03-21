#include <chrono>
#include <iomanip>

namespace DateTime::Utils
{
    std::string FormatTime(std::chrono::system_clock::time_point tp)
    {
        std::stringstream ss;
        auto t = std::chrono::system_clock::to_time_t(tp);
        auto tp2 = std::chrono::system_clock::from_time_t(t);
        if (tp2 > tp)
            t = std::chrono::system_clock::to_time_t(tp - std::chrono::seconds(1));
        ss << std::put_time(std::localtime(&t), "%Y-%m-%d %T")
           << "." << std::setfill('0') << std::setw(3)
           << (std::chrono::duration_cast<std::chrono::milliseconds>(
                   tp.time_since_epoch())
                   .count() %
               1000);
        return ss.str();
    }

    std::string CurrentTimeStr()
    {
        return FormatTime(std::chrono::system_clock::now());
    }
}
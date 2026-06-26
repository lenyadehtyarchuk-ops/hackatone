#pragma once
#include <string>
#include <vector>
#include <optional>
#include <cstdint>

struct NmeaFix {
    double timestamp_s;   // секунды от начала суток (UTC)
    double radio_alt_m;   // высота над рельефом (AGL), м
    double baro_alt_m;    // барометрическая высота (MSL), м
    int    gps_quality;   // поле 6 GPGGA: 0=нет GPS (зона глушения), 1=GPS, 2=DGPS
};

class NmeaParser {
public:
    // Разобрать одну строку GPGGA, вернуть Fix если успешно
    static std::optional<NmeaFix> parse_gpgga(const std::string& line);

    // Прочитать весь файл NMEA, вернуть все успешно разобранные фиксы
    static std::vector<NmeaFix> load_file(const std::string& path);

private:
    static bool verify_checksum(const std::string& line);
    static double parse_hhmmss(const std::string& token);
};

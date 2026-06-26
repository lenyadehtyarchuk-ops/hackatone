#include "nmea_parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <iostream>

static std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> fields;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        fields.push_back(tok);
    return fields;
}

bool NmeaParser::verify_checksum(const std::string& line) {
    // Формат: $<data>*HH
    if (line.size() < 4 || line[0] != '$') return false;
    auto star = line.rfind('*');
    if (star == std::string::npos || star + 3 > line.size()) return false;

    uint8_t cs = 0;
    for (size_t i = 1; i < star; ++i)
        cs ^= static_cast<uint8_t>(line[i]);

    unsigned int expected = 0;
    if (sscanf(line.c_str() + star + 1, "%2X", &expected) != 1) return false;
    return cs == static_cast<uint8_t>(expected);
}

double NmeaParser::parse_hhmmss(const std::string& tok) {
    if (tok.size() < 6) return 0.0;
    int hh = std::stoi(tok.substr(0, 2));
    int mm = std::stoi(tok.substr(2, 2));
    double ss = std::stod(tok.substr(4));
    return hh * 3600.0 + mm * 60.0 + ss;
}

std::optional<NmeaFix> NmeaParser::parse_gpgga(const std::string& raw_line) {
    std::string line = raw_line;
    // Обрезать \r\n
    while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
        line.pop_back();

    if (line.empty() || line[0] != '$') return std::nullopt;
    if (!verify_checksum(line)) return std::nullopt;

    // Убрать $, *XX
    auto star = line.rfind('*');
    std::string body = line.substr(1, star - 1);
    auto fields = split_csv(body);

    // GPGGA fields:
    // 0: GPGGA
    // 1: hhmmss.sss
    // 2: lat, 3: N/S, 4: lon, 5: E/W
    // 6: fix quality, 7: satellites
    // 8: HDOP
    // 9: altitude (antenna), 10: M
    // 11: geoid sep (используем как baro_alt), 12: M
    // 13: DGPS age, 14: DGPS ID

    if (fields.size() < 13) return std::nullopt;
    if (fields[0] != "GPGGA" && fields[0] != "GNGGA") return std::nullopt;

    NmeaFix fix{};
    try {
        fix.timestamp_s = fields[1].empty() ? 0.0 : parse_hhmmss(fields[1]);
        fix.gps_quality = fields[6].empty() ? 1   : std::stoi(fields[6]);
        fix.radio_alt_m = fields[9].empty() ? 0.0 : std::stod(fields[9]);
        fix.baro_alt_m  = fields[11].empty() ? 0.0 : std::stod(fields[11]);
        // Координаты присутствуют только когда GPS работает (quality > 0)
        fix.lat = 0.0; fix.lon = 0.0;
        if (fix.gps_quality > 0 && !fields[2].empty() && !fields[4].empty()) {
            // NMEA: DDMM.MMMMM → градусы
            auto nmea2deg = [](const std::string& s) {
                double raw = std::stod(s);
                int deg = static_cast<int>(raw / 100);
                return deg + (raw - deg * 100.0) / 60.0;
            };
            fix.lat = nmea2deg(fields[2]) * (fields[3] == "S" ? -1 : 1);
            fix.lon = nmea2deg(fields[4]) * (fields[5] == "W" ? -1 : 1);
        }
    } catch (...) {
        return std::nullopt;
    }

    return fix;
}

std::vector<NmeaFix> NmeaParser::load_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open NMEA file: " + path);

    std::vector<NmeaFix> result;
    std::string line;
    int line_num = 0, parsed = 0;
    while (std::getline(f, line)) {
        ++line_num;
        auto fix = parse_gpgga(line);
        if (fix) {
            result.push_back(*fix);
            ++parsed;
        }
    }
    std::cerr << "[NMEA] Прочитано строк: " << line_num
              << ", распознано: " << parsed << "\n";
    return result;
}

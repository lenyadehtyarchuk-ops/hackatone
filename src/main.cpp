#include <iostream>
#include <vector>
#include "gdal_priv.h"

int main() {
    GDALAllRegister();

    GDALDataset* poDataset = (GDALDataset*) GDALOpen("input_data.tif", GA_ReadOnly);
    if (poDataset == nullptr) {
        std::cerr << "Не удалось открыть файл карты!" << std::endl;
        return 1;
    }

    int width = poDataset->GetRasterXSize();
    int height = poDataset->GetRasterYSize();
    std::cout << "Размер сетки карты: " << width << " x " << height << " точек." << std::endl;

    GDALRasterBand* poBand = poDataset->GetRasterBand(1);

    std::vector<float> elevationMatrix(width * height);

    poBand->RasterIO(GF_Read, 0, 0, width, height,
                     &elevationMatrix[0], width, height, GDT_Float32, 0, 0);

    int centerX = width / 2;
    int centerY = height / 2;

    float centerElevation = elevationMatrix[centerY * width + centerX];

    std::cout << "Высота в центре карты: " << centerElevation << " метров над уровнем моря." << std::endl;

    GDALClose(poDataset);
    return 0;
}

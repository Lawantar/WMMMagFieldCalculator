#include <stdio.h>
#include "geomagcalc.h"

int main() {
    float h = 100000.0f;
    float lat = 0.0f;
    float lon = 120.0f;
    float dyear = ConvertSecsToDecimalYear(773176800); // 02.07.2024 в секундах от 01.01.2000
    Vector test = GetMagFieldElements(dyear, lat, lon, h);
    printf("%f\n", dyear);
    printf("North Comp - X (nT): %f\n", test.x);
    printf("East Comp - Y (nT): %f\n", test.y);
    printf("Vertical Comp - Z (nT): %f\n", test.z);
}

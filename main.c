#include <stdio.h>
#include "geomagcalc.h"

int main() {
    float h = 100000.0f;
    float lat = 0.0f;
    float lon = 120.0f;
    float dyear = 2022.5f;
    Elements test = GetMagFieldElements(dyear, lat, lon, h);
    printf("Declination (deg): %f\n", test.declination);
    printf("Inclination (deg): %f\n", test.inclination);
    printf("Horizontal Intensity (nT): %f\n", test.horizontal);
    printf("North Comp - X (nT): %f\n", test.north);
    printf("East Comp - Y (nT): %f\n", test.east);
    printf("Vertical Comp - Z (nT): %f\n", test.vertical);
    printf("Total Field (nT): %f\n", test.total);
}

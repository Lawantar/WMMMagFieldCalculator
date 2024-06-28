#include <stdio.h>
#include "geomagcalc.h"

int main() {
    float h = 100000.0f;
    float lat = 55.031667f;
    float lon = 82.927778f;
    float dyear = 2024.49f;
    Vector position = ConvertGeodeticToEcef(lat, lon, h);
    Vector mag_field = GeoMag(dyear, position);
    Elements out = ConvertMagFieldToElements(mag_field, lat, lon);
    printf("Declination (deg): %f\n", out.declination);
    printf("Inclination (deg): %f\n", out.inclination);
    printf("Horizontal Intensity (nT): %f\n", out.horizontal);
    printf("North Comp - X (nT): %f\n", out.north);
    printf("East Comp - Y (nT): %f\n", out.east);
    printf("Vertical Comp - Z (nT): %f\n", out.vertical);
    printf("Total Field (nT): %f\n", out.total);
}

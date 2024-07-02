#include "geomagcalc.h"
#include "math.h"

F_32 GetConstModelC(SI_16 n, SI_16 m, F_32 dyear) {
    SI_16 index = (m * (2 * NMAX - m + 1)) / 2 + n;
    return WMM2020.Main_Field_Coeff_C[index] + (dyear - WMM2020.epoch) * WMM2020.Secular_Var_Coeff_C[index];
}

F_32 GetConstModelS(SI_16 n, SI_16 m, F_32 dyear) {
    SI_16 index = (m * (2 * NMAX - m + 1)) / 2 + n;
    return WMM2020.Main_Field_Coeff_S[index] + (dyear - WMM2020.epoch) * WMM2020.Secular_Var_Coeff_S[index];
}

// Возвращает положение в координатах ITRS (метры)
// Аргументы:
// lat - Широта в градусах
// lon - Долгота в градусах
// height - Выоста над поверхностью в системе WGS 84 ellipsoid в метрах

Vector ConvertGeodeticToEcef(F_32 lat, F_32 lon, F_32 height) {
    Vector position;
    F_32 phi = lat * ((F_32) (M_PI / 180.0));
    F_32 lam = lon * ((F_32) (M_PI / 180.0));
    const F_32 a = 6378137.0f;
    const F_32 e2 = 0.00669438f;
    const F_32 e2m = 0.99330562f;
    F_32 sphi = sinf(phi);
    F_32 cphi = cosf(phi);
    F_32 slam = sinf(lam);
    F_32 clam = cosf(lam);
    F_32 n = a / sqrtf(1.0f - e2 * (sphi * sphi));
    F_32 z = (e2m * n + height) * sphi;
    F_32 r = (n + height) * cphi;
    position.x = r * clam;
    position.y = r * slam;
    position.z = z;
    return position;
}

// Возвращает вектор магнитного поля в ITRS (нТл)
// Аргументы:
// dyear - Дата в виде десятичной дроби (2020.0 - 1 января 2020 или 2022.5 - 1 июля 2022)
// position - Положение в координатах ITRS

Vector GeoMag(F_32 dyear, Vector position) {
    Vector GeoMag;
    F_32 x = position.x;
    F_32 y = position.y;
    F_32 z = position.z;
    F_32 px = 0.0f;
    F_32 py = 0.0f;
    F_32 pz = 0.0f;
    F_32 rsqrd = x * x + y * y + z * z;
    F_32 temp = EARTH_R / rsqrd;
    F_32 a = x * temp;
    F_32 b = y * temp;
    F_32 f = z * temp;
    F_32 g = EARTH_R * temp;

    SI_16 n, m;
    F_32 Vtop = EARTH_R / sqrtf(rsqrd);
    F_32 Wtop = 0.0f;
    F_32 Vprev = 0.0f;
    F_32 Wprev = 0.0f;
    F_32 Vnm = Vtop;
    F_32 Wnm = Wtop;

    for (m = 0; m <= NMAX + 1; ++m) {
        for (n = m; n <= NMAX + 1; ++n) {
            if (n == m) {
                if (m != 0) {
                    temp = Vtop;
                    Vtop = (F_32)(2 * m - 1) * (a * Vtop - b * Wtop);
                    Wtop = (F_32)(2 * m - 1) * (a * Wtop + b * temp);
                    Vprev = 0;
                    Wprev = 0;
                    Vnm = Vtop;
                    Wnm = Wtop;
                }
            } else {
                temp = Vnm;
                F_32 invs_temp = 1.0f / ((F_32) (n - m));
                Vnm = ((F_32)(2 * n - 1) * f * Vnm - (F_32)(n + m - 1) * g * Vprev) * invs_temp;
                Vprev = temp;
                temp = Wnm;
                Wnm = ((F_32)(2 * n - 1) * f * Wnm - (F_32)(n + m - 1) * g * Wprev) * invs_temp;
                Wprev = temp;
            }
            if (m < NMAX && n >= m + 2) {
                px += 0.5f * (F_32)(n - m) * (F_32)(n - m - 1) * (GetConstModelC(n - 1, m + 1, dyear) * Vnm + GetConstModelS(n - 1, m + 1, dyear) * Wnm);
                py += 0.5f * (F_32)(n - m) * (F_32)(n - m - 1) * (-GetConstModelC(n - 1, m + 1, dyear) * Wnm + GetConstModelS(n - 1, m + 1, dyear) * Vnm);
            }
            if (n >= 2 && m >= 2) {
                px += 0.5f * (-GetConstModelC(n - 1, m - 1, dyear) * Vnm - GetConstModelS(n - 1, m - 1, dyear) * Wnm);
                py += 0.5f * (-GetConstModelC(n - 1, m - 1, dyear) * Wnm + GetConstModelS(n - 1, m - 1, dyear) * Vnm);
            }
            if (m == 1 && n >= 2) {
                px += -GetConstModelC(n - 1, 0, dyear) * Vnm;
                py += -GetConstModelC(n - 1, 0, dyear) * Wnm;
            }
            if (n >= 2 && n > m) {
                pz += (F_32)(n - m) * (-GetConstModelC(n - 1, m, dyear) * Vnm - GetConstModelS(n - 1, m, dyear) * Wnm);
            }
        }
    }
    GeoMag.x = -px;
    GeoMag.y = -py;
    GeoMag.z = -pz;
    return GeoMag;
}

// Возвращает элементы магнитного поля (X,Y,Z,H,F,I,D) в нТл и градусах
// Аргументы:
// dyear - Дата в виде десятичной дроби (2020.0 - 1 января 2020 или 2022.5 - 1 июля 2022)
// lat - Широта в градусах
// lon - Долгота в градусах
// height - Выоста над поверхностью в системе WGS 84 ellipsoid в метрах

Elements GetMagFieldElements(F_32 dyear, F_32 lat, F_32 lon, F_32 height) {
    Vector position = ConvertGeodeticToEcef(lat, lon, height);
    Vector mag_field = GeoMag(dyear, position);
    Elements MagFieldElements;
    F_32 phi = lat * ((F_32) (M_PI / 180.0f));
    F_32 lam = lon * ((F_32) (M_PI / 180.0f));
    F_32 x1 = cosf(lam) * mag_field.x + sinf(lam) * mag_field.y;
    F_32 north = -sinf(phi) * x1 + cosf(phi) * mag_field.z;
    F_32 east = -sinf(lam) * mag_field.x + cosf(lam) * mag_field.y;
    F_32 vertical = -cosf(phi) * x1 + -sinf(phi) * mag_field.z;
    F_32 horizontal = sqrtf(north * north + east * east);
    MagFieldElements.north = north;
    MagFieldElements.east = east;
    MagFieldElements.vertical = vertical;
    MagFieldElements.horizontal = horizontal;
    MagFieldElements.total = sqrtf(horizontal * horizontal + vertical * vertical);
    MagFieldElements.inclination = atan2f(vertical, horizontal) * ((F_32) (180.0f / M_PI));
    MagFieldElements.declination = atan2f(east, north) * ((F_32) (180.0f / M_PI));
    return MagFieldElements;
}

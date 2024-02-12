﻿#ifndef _ASTRONOMY_H_
#define _ASTRONOMY_H_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstddef>

enum solar_term_t {
    spring_equinox = 0,
    pure_brightness = 15,
    grain_rain = 30,
    beginning_of_summer = 45,
    grain_buds = 60,
    grain_ear = 75,
    summer_solstice = 90,
    minor_heat = 105,
    major_heat = 120,
    beginning_of_autumn = 135,
    end_of_heat = 150,
    white_dew = 165,
    autumn_equinox = 180,
    cold_dew = 195,
    frost_s_descent = 210,
    beginning_of_winter = 225,
    minor_snow = 240,
    major_snow = 255,
    winter_solstice = 270,
    minor_cold = 285,
    major_cold = 300,
    beginning_of_spring = 315,
    rain_water = 330,
    awakening_of_insects = 345,
};

namespace astronomy {
    typedef long double REAL;

    static constexpr REAL JD2000 = 2451545;  // 2000年前儒略日数(2000-1-1 12:00:00格林威治平时)

    struct daytime_t {
        int year;
        int month;
        int day;
        int hour;
        int minute;
        REAL second;
    };

    struct ecliptic_position_t {
        REAL longitude, latitude;
    };

    namespace detail {
        static constexpr REAL RADIAN_PER_DEGREE = M_PI / 180.0;
        static constexpr REAL DEGREE_PER_RADIAN = 180.0 / M_PI;
        static constexpr REAL PI_2 = M_PI * 2;

        struct vsop87_coefficient_t {
            REAL a, b, c;
        };

        struct nutation_coefficient_t {
            REAL a0, a1, a2, a3, a4, sin1, sin2, cos1, cos2;
        };

        struct delta_time_t {
            int y;
            REAL a0, a1, a2, a3;
        };

        struct elp2000_coefficient_t {
            REAL f, a0, a1, a2, a3, a4;
        };

        static REAL clamp_randians(REAL a) {
            while (a < 0) a += PI_2;
            while (a > PI_2) a -= PI_2;
            return a;
        }

        static REAL vsop87_periodic_terms(const vsop87_coefficient_t *c, std::size_t n, REAL t) {
            REAL v = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const auto &e = c[i];
                v += e.a * std::cos(e.b + e.c * t);
            }
            return v;
        }

        static REAL elp2000_periodic_terms(const elp2000_coefficient_t *c, std::size_t n, REAL t) {
            REAL v = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const auto &e = c[i];
                v += e.f * std::sin(e.a0 + (e.a1 + (e.a2 + (e.a3 + e.a4 * t) * t) * t) * t);
            }
            return v;
        }

        template <class Dummy>
        struct impl {
            static const vsop87_coefficient_t E10[];
            static const vsop87_coefficient_t E11[];
            static const vsop87_coefficient_t E12[];
            static const vsop87_coefficient_t E13[];
            static const vsop87_coefficient_t E14[];
            static const vsop87_coefficient_t E15[];

            static const vsop87_coefficient_t E20[];
            static const vsop87_coefficient_t E21[];

            static const vsop87_coefficient_t E30[];
            static const vsop87_coefficient_t E31[];
            static const vsop87_coefficient_t E32[];
            static const vsop87_coefficient_t E33[];

            static const elp2000_coefficient_t M10[];
            static const elp2000_coefficient_t M11[];
            static const elp2000_coefficient_t M12[];

            static const elp2000_coefficient_t M20[];
            static const elp2000_coefficient_t M21[];

            static const nutation_coefficient_t NT[];

            static const delta_time_t D[];

            static REAL calc_earth_longitude(REAL t);
            static REAL calc_earth_latitude(REAL t);

            static void adjust_sun_aberration_and_nutation(REAL t, ecliptic_position_t &pos);
            static void adjust_sun_aberration_and_nutation_2(REAL t, ecliptic_position_t &pos);
            static ecliptic_position_t calc_sun_position(REAL jd);

            static REAL calc_moon_longitude(REAL t);
            static REAL calc_moon_latitude(REAL t);

            static REAL adjust_precession(REAL jd, REAL l);
            static REAL calc_moon_ecliptic_longitude(REAL jd);
            static REAL get_moon_ecliptic_longitude(REAL jd);

            static REAL get_sun_ecliptic_longitude(REAL jd);

            static REAL calc_delta_t(REAL t);
        };

        // 地球运动VSOP87参数 黄经周期项
        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E10[] = {
            { 175347045673, 0            ,      0            },
            {   3341656456, 4.66925680417,   6283.0758499914 },
            {     34894275, 4.62610241759,  12566.1516999828 },
            {      3417571, 2.82886579606,      3.5231183490 },
            {      3497056, 2.74411800971,   5753.3848848968 },
            {      3135896, 3.62767041758,  77713.7714681205 },
            {      2676218, 4.41808351397,   7860.4193924392 },
            {      2342687, 6.13516237631,   3930.2096962196 },
            {      1273166, 2.03709655772,    529.6909650946 },
            {      1324292, 0.74246356352,  11506.7697697936 },
            {       901855, 2.04505443513,     26.2983197998 },
            {      1199167, 1.10962944315,   1577.3435424478 },
            {       857223, 3.50849156957,    398.1490034082 },
            {       779786, 1.17882652114,   5223.6939198022 },
            {       990250, 5.23268129594,   5884.9268465832 },
            {       753141, 2.53339053818,   5507.5532386674 },
            {       505264, 4.58292563052,  18849.2275499742 },
            {       492379, 4.20506639861,    775.5226113240 },
            {       356655, 2.91954116867,      0.0673103028 },
            {       284125, 1.89869034186,    796.2980068164 },
            {       242810, 0.34481140906,   5486.7778431750 },
            {       317087, 5.84901952218,  11790.6290886588 },
            {       271039, 0.31488607649,  10977.0788046990 },
            {       206160, 4.80646606059,   2544.3144198834 },
            {       205385, 1.86947813692,   5573.1428014331 },
            {       202261, 2.45767795458,   6069.7767545534 },
            {       126184, 1.08302630210,     20.7753954924 },
            {       155516, 0.83306073807,    213.2990954380 },
            {       115132, 0.64544911683,      0.9803210682 },
            {       102851, 0.63599846727,   4694.0029547076 },
            {       101724, 4.26679821365,      7.1135470008 },
            {        99206, 6.20992940258,   2146.1654164752 },
            {       132212, 3.41118275555,   2942.4634232916 },
            {        97607, 0.68101272270,    155.4203994342 },
            {        85128, 1.29870743025,   6275.9623029906 },
            {        74651, 1.75508916159,   5088.6288397668 },
            {       101895, 0.97569221824,  15720.8387848784 },
            {        84711, 3.67080093025,  71430.6956181291 },
            {        73547, 4.67926565481,    801.8209311238 },
            {        73874, 3.50319443167,   3154.6870848956 },
            {        78756, 3.03698313141,  12036.4607348882 },
            {        79637, 1.80791330700,  17260.1546546904 },
            {        85803, 5.98322631256, 161000.6857376741 },
            {        56963, 2.78430398043,   6286.5989683404 },
            {        61148, 1.81839811024,   7084.8967811152 },
            {        69627, 0.83297596966,   9437.7629348870 },
            {        56116, 4.38694880779,  14143.4952424306 },
            {        62449, 3.97763880587,   8827.3902698748 },
            {        51145, 0.28306864501,   5856.4776591154 },
            {        55577, 3.47006009062,   6279.5527316424 },
            {        41036, 5.36817351402,   8429.2412664666 },
            {        51605, 1.33282746983,   1748.0164130670 },
            {        51992, 0.18914945834,  12139.5535091068 },
            {        49000, 0.48735065033,   1194.4470102246 },
            {        39200, 6.16832995016,  10447.3878396044 },
            {        35566, 1.77597314691,   6812.7668150860 },
            {        36770, 6.04133859347,  10213.2855462110 },
            {        36596, 2.56955238628,   1059.3819301892 },
            {        33291, 0.59309499459,  17789.8456197850 },
            {        35954, 1.70876111898,   2352.8661537718 },
        };

        // 黄经泊松1项
        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E11[] = {
            { 628331966747491, 0            ,     0            },
            {       206058863, 2.67823455584,  6283.0758499914 },
            {         4303430, 2.63512650414, 12566.1516999828 },
            {          425264, 1.59046980729,     3.5231183490 },
            {          108977, 2.96618001993,  1577.3435424478 },
            {           93478, 2.59212835365, 18849.2275499742 },
            {          119261, 5.79557487799,    26.2983197998 },
            {           72122, 1.13846158196,   529.6909650946 },
            {           67768, 1.87472304791,   398.1490034082 },
            {           67327, 4.40918235168,  5507.5532386674 },
            {           59027, 2.88797038460,  5223.6939198022 },
            {           55976, 2.17471680261,   155.4203994342 },
            {           45407, 0.39803079805,   796.2980068164 },
            {           36369, 0.46624739835,   775.5226113240 },
            {           28958, 2.64707383882,     7.1135470008 },
            {           19097, 1.84628332577,  5486.7778431750 },
            {           20844, 5.34138275149,     0.9803210682 },
            {           18508, 4.96855124577,   213.2990954380 },
            {           16233, 0.03216483047,  2544.3144198834 },
            {           17293, 2.99116864949,  6275.9623029906 },
        };

        // 黄经泊松2项
        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E12[] = {
            { 52918870, 0            ,     0            },
            {  8719837, 1.07209665242,  6283.0758499914 },
            {   309125, 0.86728818832, 12566.1516999828 },
            {    27339, 0.05297871691,     3.5231183490 },
            {    16334, 5.18826691036,    26.2983197998 },
            {    15752, 3.68457889430,   155.4203994342 },
            {     9541, 0.75742297675, 18849.2275499742 },
            {     8937, 2.05705419118, 77713.7714681205 },
            {     6952, 0.82673305410,   775.5226113240 },
            {     5064, 4.66284525271,  1577.3435424478 }
        };

        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E13[] = {
            { 289226, 5.84384198723,  6283.0758499914 },
            {  34955, 0            ,     0            },
            {  16819, 5.48766912348, 12566.1516999828 },
        };

        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E14[] = {
            { 114084, 3.14159265359,     0            },
            {   7717, 4.13446589358,  6283.0758499914 },
            {    765, 3.83803776214, 12566.1516999828 },
        };

        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E15[] = {
            { 878, 3.14159265359, 0 },
        };

        // 日心黄经
        template <class Dummy>
        REAL impl<Dummy>::calc_earth_longitude(REAL t) {
            REAL L0 = vsop87_periodic_terms(E10, sizeof(E10) / sizeof(*E10), t);
            REAL L1 = vsop87_periodic_terms(E11, sizeof(E11) / sizeof(*E11), t);
            REAL L2 = vsop87_periodic_terms(E12, sizeof(E12) / sizeof(*E12), t);
            REAL L3 = vsop87_periodic_terms(E13, sizeof(E13) / sizeof(*E13), t);
            REAL L4 = vsop87_periodic_terms(E14, sizeof(E14) / sizeof(*E14), t);
            REAL L5 = vsop87_periodic_terms(E15, sizeof(E15) / sizeof(*E15), t);

            return (L0 + (L1 + (L2 + (L3 + (L4 + L5 * t) * t) * t) * t) * t) / 1E11;
        }

        // 黄纬周期项
        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E20[] = {
            { 279620, 3.19870156017, 84334.6615813083 },
            { 101643, 5.42248619256,  5507.5532386674 },
            {  80445, 3.88013204458,  5223.6939198022 },
            {  43806, 3.70444689758,  2352.8661537718 },
            {  31933, 4.00026369781,  1577.3435424478 },
            {  22724, 3.98473831560,  1047.7473117547 },
            {  16392, 3.56456119782,  5856.4776591154 },
            {  18141, 4.98367470263,  6283.0758499914 },
            {  14443, 3.70275614914,  9437.7629348870 },
            {  14304, 3.41117857525, 10213.2855462110 },
        };

        template <class Dummy>
        const vsop87_coefficient_t impl<Dummy>::E21[] = {
            { 9030, 3.89729061890, 5507.5532386674 },
            { 6177, 1.73038850355, 5223.6939198022 },
        };

        // 日心黄纬
        template <class Dummy>
        REAL impl<Dummy>::calc_earth_latitude(REAL t) {
            REAL L0 = vsop87_periodic_terms(E20, sizeof(E20) / sizeof(*E20), t);
            REAL L1 = vsop87_periodic_terms(E21, sizeof(E21) / sizeof(*E21), t);

            return (L0 + L1 * t) / 1E11;
        }

        template <class Dummy>
        const nutation_coefficient_t impl<Dummy>::NT[] = {
            { 2.1824391966,   -33.757045954,  0.0000362262, 3.7340E-08, -2.8793E-10, -171996, -1742, 92025,  89 },
            { 3.5069406862,  1256.663930738,  0.0000105845, 6.9813E-10, -2.2815E-10,  -13187,   -16,  5736, -31 },
            { 1.3375032491, 16799.418221925, -0.0000511866, 6.4626E-08, -5.3543E-10,   -2274,    -2,   977,  -5 },
            { 4.3648783932,   -67.514091907,  0.0000724525, 7.4681E-08, -5.7586E-10,    2062,     2,  -895,   5 },
            { 0.0431251803,  -628.301955171,  0.0000026820, 6.5935E-10,  5.5705E-11,   -1426,    34,    54,  -1 },
            { 2.3555557435,  8328.691425719,  0.0001545547, 2.5033E-07, -1.1863E-09,     712,     1,    -7,   0 },
            { 3.4638155059,  1884.965885909,  0.0000079025, 3.8785E-11, -2.8386E-10,    -517,    12,   224,  -6 },
            { 5.4382493597, 16833.175267879, -0.0000874129, 2.7285E-08, -2.4750E-10,    -386,    -4,   200,   0 },
            { 3.6930589926, 25128.109647645,  0.0001033681, 3.1496E-07, -1.7218E-09,    -301,     0,   129,  -1 },
            { 3.5500658664,   628.361975567,  0.0000132664, 1.3575E-09, -1.7245E-10,     217,    -5,   -95,   3 },
        };

        // 光行差与天体章动
        // 注释掉的部分为黄纬，计算农历不需要黄纬
        template <class Dummy>
        void impl<Dummy>::adjust_sun_aberration_and_nutation(REAL t, ecliptic_position_t &pos) {
            static const REAL E[] = { 0.016708634, -0.000042037, -0.0000001267 };  // 离心率
            static const REAL P[] = { 102.93735 * RADIAN_PER_DEGREE, 1.71946 * RADIAN_PER_DEGREE, 0.00046 * RADIAN_PER_DEGREE };  // 近点
            static const REAL L[] = { 280.4664567 * RADIAN_PER_DEGREE, 36000.76982779 * RADIAN_PER_DEGREE, 0.0003032028 * RADIAN_PER_DEGREE, RADIAN_PER_DEGREE / 49931000.0, RADIAN_PER_DEGREE / -153000000.0 };  // 太平黄经
            static const REAL K = 20.49552 * RADIAN_PER_DEGREE / 3600.0;  // 光行差常数，单位角秒，这里转化为弧度

            // 光行差
            REAL t1 = t / 36525.0;
            REAL l = L[0] + (L[1] + (L[2] + (L[3] + L[4] * t1) * t1) * t1) * t1;
            REAL p = P[0] + (P[1] + P[2] * t1) * t1;
            REAL e = E[0] + (E[1] + E[2] * t1) * t1;
            REAL dL = l - pos.longitude;
            REAL dP = p - pos.longitude;

            pos.longitude -= K * (std::cos(dL) - e * std::cos(dP)) / std::cos(pos.latitude);
            //pos.latitude -= K * std::sin(pos.latitude) * (std::sin(dL) - e * std::sin(dP));

            pos.longitude = clamp_randians(pos.longitude);

            // 天体章动
            REAL longitude = 0;
            //REAL obliquity = 0;
            for (std::size_t i = 0, c = sizeof(NT) / sizeof(*NT); i < c; ++i) {
                const auto &n = NT[i];
                REAL v = n.a0 + (n.a1 + (n.a2 + (n.a3 + n.a4 * t1) * t1) * t1) * t1;
                longitude += (n.sin1 + n.sin2 * t1 / 10) * std::sin(v);
                //obliquity += (n.cos1 + n.cos2 * t1 / 10) * std::cos(v);
            }

            longitude /= (36000000.0 * DEGREE_PER_RADIAN);
            //obliquity /= (36000000.0 * DEGREE_PER_RADIAN);
            pos.longitude = clamp_randians(pos.longitude + longitude);
        }

        // 太阳视位置
        template <class Dummy>
        ecliptic_position_t impl<Dummy>::calc_sun_position(REAL jd) {
            REAL t = jd / 365250.0;

            ecliptic_position_t pos{};
            pos.longitude = clamp_randians(calc_earth_longitude(t) + M_PI);  // 地心黄经 = 日心黄经 + 180度
            pos.latitude = -calc_earth_latitude(t);  // 地心黄纬 = -日心黄纬
            adjust_sun_aberration_and_nutation(jd, pos);  // 修正天体章动
            return pos;
        }

        template <class Dummy>
        const elp2000_coefficient_t impl<Dummy>::M10[] = {
            { 22639.5858800,   2.3555545723,   8328.6914247251,  1.5231275E-04,  2.5041111E-07, -1.1863391E-09 },
            {  4586.4383203,   8.0413790709,   7214.0628654588, -2.1850087E-04, -1.8646419E-07,  8.7760973E-10 },
            {  2369.9139357,  10.3969336431,  15542.7542901840, -6.6188121E-05,  6.3946925E-08, -3.0872935E-10 },
            {   769.0257187,   4.7111091445,  16657.3828494503,  3.0462550E-04,  5.0082223E-07, -2.3726782E-09 },
            {  -666.4175399,  -0.0431256817,    628.3019552485, -2.6638815E-06,  6.1639211E-10, -5.4439728E-11 },
            {  -411.5957339,   3.2558104895,  16866.9323152810, -1.2804259E-04, -9.8998954E-09,  4.0433461E-11 },
            {   211.6555524,   5.6858244986,  -1114.6285592663, -3.7081362E-04, -4.3687530E-07,  2.0639488E-09 },
            {   205.4359530,   8.0845047526,   6585.7609102104, -2.1583699E-04, -1.8708058E-07,  9.3204945E-10 },
            {   191.9561973,  12.7524882154,  23871.4457149091,  8.6124629E-05,  3.1435804E-07, -1.4950684E-09 },
            {   164.7286185,  10.4400593249,  14914.4523349355, -6.3524240E-05,  6.3330532E-08, -2.5428962E-10 },
            {  -147.3213842,  -2.3986802540,  -7700.3894694766, -1.5497663E-04, -2.4979472E-07,  1.1318993E-09 },
            {  -124.9881185,   5.1984668216,   7771.3771450920, -3.3094061E-05,  3.1973462E-08, -1.5436468E-10 },
            {  -109.3803637,   2.3124288905,   8956.9933799736,  1.4964887E-04,  2.5102751E-07, -1.2407788E-09 },
            {    55.1770578,   7.1411231536,  -1324.1780250970,  6.1854469E-05,  7.3846820E-08, -3.4916281E-10 },
            {   -45.0996092,   5.6113650618,  25195.6237400061,  2.4270161E-05,  2.4051122E-07, -1.1459056E-09 },
            {    39.5333010,  -0.9002559173,  -8538.2408905558,  2.8035534E-04,  2.6031101E-07, -1.2267725E-09 },
            {    38.4298346,  18.4383127140,  22756.8171556428, -2.8468899E-04, -1.2251727E-07,  5.6888037E-10 },
            {    36.1238141,   7.0666637168,  24986.0742741754,  4.5693825E-04,  7.5123334E-07, -3.5590172E-09 },
            {    30.7725751,  16.0827581417,  14428.1257309177, -4.3700174E-04, -3.7292838E-07,  1.7552195E-09 },
            {   -28.3971008,   7.9982533891,   7842.3648207073, -2.2116475E-04, -1.8584780E-07,  8.2317000E-10 },
            {   -24.3582283,  10.3538079614,  16171.0562454324, -6.8852003E-05,  6.4563317E-08, -3.6316908E-10 },
            {   -18.5847068,   2.8429122493,   -557.3142796331, -1.8540681E-04, -2.1843765E-07,  1.0319744E-09 },
            {    17.9544674,   5.1553411398,   8399.6791003405, -3.5757942E-05,  3.2589854E-08, -2.0880440E-10 },
            {    14.5302779,  12.7956138971,  23243.1437596606,  8.8788511E-05,  3.1374165E-07, -1.4406287E-09 },
            {    14.3796974,  15.1080427876,  32200.1371396342,  2.3843738E-04,  5.6476915E-07, -2.6814075E-09 },
            {    14.2514576, -24.0810366320,     -2.3011998397,  1.5231275E-04,  2.5041111E-07, -1.1863391E-09 },
            {    13.8990596,  20.7938672862,  31085.5085803679, -1.3237624E-04,  1.2789385E-07, -6.1745870E-10 },
            {    13.1940636,   3.3302699264,  -9443.3199839914, -5.2312637E-04, -6.8728642E-07,  3.2502879E-09 },
            {    -9.6790568,  -4.7542348263, -16029.0808942018, -3.0728938E-04, -5.0020584E-07,  2.3182384E-09 },
            {    -9.3658635,  11.2971895604,  24080.9951807398, -3.4654346E-04, -1.9636409E-07,  9.1804319E-10 },
            {     8.6055318,   5.7289501804,  -1742.9305145148, -3.6814974E-04, -4.3749170E-07,  2.1183885E-09 },
            {    -8.4530982,   7.5540213938,  16100.0685698171,  1.1921869E-04,  2.8238458E-07, -1.3407038E-09 },
            {     8.0501724,  10.4831850066,  14286.1503796870, -6.0860358E-05,  6.2714140E-08, -1.9984990E-10 },
            {    -7.6301553,   4.6679834628,  17285.6848046987,  3.0196162E-04,  5.0143862E-07, -2.4271179E-09 },
            {    -7.4474952,  -0.0862513635,   1256.6039104970, -5.3277630E-06,  1.2327842E-09, -1.0887946E-10 },
            {     7.3712011,   8.1276304344,   5957.4589549619, -2.1317311E-04, -1.8769697E-07,  9.8648918E-10 },
            {     7.0629900,   0.9591375719,     33.7570471374, -3.0829302E-05, -3.6967043E-08,  1.7385419E-10 },
            {    -6.3831491,   9.4966777258,   7004.5133996281,  2.1416722E-04,  3.2425793E-07, -1.5355019E-09 },
            {    -5.7416071,  13.6527441326,  32409.6866054649, -1.9423071E-04,  5.4047029E-08, -2.6829589E-10 },
            {     4.3740095,  18.4814383957,  22128.5152003943, -2.8202511E-04, -1.2313366E-07,  6.2332010E-10 },
            {    -3.9976134,   7.9669196340,  33524.3151647312,  1.7658291E-04,  4.9092233E-07, -2.3322447E-09 },
            {    -3.2096876,  13.2398458924,  14985.4400105508, -2.5159493E-04, -1.5449073E-07,  7.2324505E-10 },
            {    -2.9145404,  12.7093625336,  24499.7476701576,  8.3460748E-05,  3.1497443E-07, -1.5495082E-09 },
            {     2.7318890,  16.1258838235,  13799.8237756692, -4.3433786E-04, -3.7354477E-07,  1.8096592E-09 },
            {    -2.5679459,  -2.4418059357,  -7072.0875142282, -1.5764051E-04, -2.4917833E-07,  1.0774596E-09 },
            {    -2.5211990,   7.9551277074,   8470.6667759558, -2.2382863E-04, -1.8523141E-07,  7.6873027E-10 },
            {     2.4888871,   5.6426988169,   -486.3266040178, -3.7347750E-04, -4.3625891E-07,  2.0095091E-09 },
            {     2.1460741,   7.1842488353,  -1952.4799803455,  6.4518350E-05,  7.3230428E-08, -2.9472308E-10 },
            {     1.9777270,  23.1494218585,  39414.2000050930,  1.9936508E-05,  3.7830496E-07, -1.8037978E-09 },
            {     1.9336825,   9.4222182890,  33314.7656989005,  6.0925100E-04,  1.0016445E-06, -4.7453563E-09 },
            {     1.8707647,  20.8369929680,  30457.2066251194, -1.2971236E-04,  1.2727746E-07, -5.6301898E-10 },
            {    -1.7529659,   0.4873576771,  -8886.0057043583, -3.3771956E-04, -4.6884877E-07,  2.2183135E-09 },
            {    -1.4371624,   7.0979974718,   -695.8760698485,  5.9190587E-05,  7.4463212E-08, -4.0360254E-10 },
            {    -1.3725701,   1.4552986550,   -209.5494658307,  4.3266809E-04,  5.1072212E-07, -2.4131116E-09 },
            {     1.2618162,   7.5108957121,  16728.3705250656,  1.1655481E-04,  2.8300097E-07, -1.3951435E-09 },
        };

        template <class Dummy>
        const elp2000_coefficient_t impl<Dummy>::M11[] = {
            { 1.6768000,  -0.0431256817,   628.3019552485, -2.6638815E-06,  6.1639211E-10, -5.4439728E-11 },
            { 0.5164200,  11.2260974062,  6585.7609102104, -2.1583699E-04, -1.8708058E-07,  9.3204945E-10 },
            { 0.4138300,  13.5816519784, 14914.4523349355, -6.3524240E-05,  6.3330532E-08, -2.5428962E-10 },
            { 0.3711500,   5.5402729076,  7700.3894694766,  1.5497663E-04,  2.4979472E-07, -1.1318993E-09 },
            { 0.2756000,   2.3124288905,  8956.9933799736,  1.4964887E-04,  2.5102751E-07, -1.2407788E-09 },
            { 0.2459863, -25.6198212459,    -2.3011998397,  1.5231275E-04,  2.5041111E-07, -1.1863391E-09 },
            { 0.0711800,   7.9982533891,  7842.3648207073, -2.2116475E-04, -1.8584780E-07,  8.2317000E-10 },
            { 0.0612800,  10.3538079614, 16171.0562454324, -6.8852003E-05,  6.4563317E-08, -3.6316908E-10 },
        };

        template <class Dummy>
        const elp2000_coefficient_t impl<Dummy>::M12[] = {
            { 0.0048700,  -0.0431256817,  628.3019552485, -2.6638815E-06,  6.1639211E-10, -5.4439728E-11 },
            { 0.0022800, -27.1705318325,   -2.3011998397,  1.5231275E-04,  2.5041111E-07, -1.1863391E-09 },
            { 0.0015000,  11.2260974062, 6585.7609102104, -2.1583699E-04, -1.8708058E-07,  9.3204945E-10 },
        };

        template <class Dummy>
        REAL impl<Dummy>::calc_moon_longitude(REAL t) {
            //月球平黄经系数
            static const REAL E[] = { 3.81034392032, 8.39968473021E+03, -3.31919929753E-05, 3.20170955005E-08, -1.53637455544E-10 };

            // 岁差
            //static const REAL P[] = { 0, 50287.92262, 111.24406, 0.07699, -0.23479, -0.00178, 0.00018, 0.00001 };

            REAL L0 = elp2000_periodic_terms(M10, sizeof(M10) / sizeof(*M10), t);
            REAL L1 = elp2000_periodic_terms(M11, sizeof(M11) / sizeof(*M11), t);
            REAL L2 = elp2000_periodic_terms(M12, sizeof(M12) / sizeof(*M12), t);

            REAL L = L0 + (L1 + L2 * t) * t;
            L *= (RADIAN_PER_DEGREE / 3600);
            L += E[0] + (E[1] + (E[2] + (E[3] + E[4] * t) * t) * t) * t;

            return clamp_randians(L);
        }

        template <class Dummy>
        const elp2000_coefficient_t impl<Dummy>::M20[] = {
            18461.2400600,  1.6279052448,   8433.4661576405, -6.4021295E-05, -4.9499477E-09,  2.0216731E-11,
             1010.1671484,  3.9834598170,  16762.1575823656,  8.8291456E-05,  2.4546117E-07, -1.1661223E-09,
              999.6936555,  0.7276493275,   -104.7747329154,  2.1633405E-04,  2.5536106E-07, -1.2065558E-09,
              623.6524746,  8.7690283983,   7109.2881325435, -2.1668263E-06,  6.8896872E-08, -3.2894608E-10,
              199.4837596,  9.6692843156,  15647.5290230993, -2.8252217E-04, -1.9141414E-07,  8.9782646E-10,
              166.5741153,  6.4134738261,  -1219.4032921817, -1.5447958E-04, -1.8151424E-07,  8.5739300E-10,
              117.2606951, 12.0248388879,  23976.2204478244, -1.3020942E-04,  5.8996977E-08, -2.8851262E-10,
               61.9119504,  6.3390143893,  25090.8490070907,  2.4060421E-04,  4.9587228E-07, -2.3524614E-09,
               33.3572027, 11.1245829706,  15437.9795572686,  1.5014592E-04,  3.1930799E-07, -1.5152852E-09,
               31.7596709,  3.0832038997,   8223.9166918098,  3.6864680E-04,  5.0577218E-07, -2.3928949E-09,
               29.5766003,  8.8121540801,   6480.9861772950,  4.9705523E-07,  6.8280480E-08, -2.7450635E-10,
               15.5662654,  4.0579192538,  -9548.0947169068, -3.0679233E-04, -4.3192536E-07,  2.0437321E-09,
               15.1215543, 14.3803934601,  32304.9118725496,  2.2103334E-05,  3.0940809E-07, -1.4748517E-09,
              -12.0941511,  8.7259027166,   7737.5900877920, -4.8307078E-06,  6.9513264E-08, -3.8338581E-10,
                8.8681426,  9.7124099974,  15019.2270678508, -2.7985829E-04, -1.9203053E-07,  9.5226618E-10,
                8.0450400,  0.6687636586,   8399.7091105030, -3.3191993E-05,  3.2017096E-08, -1.5363746E-10,
                7.9585542, 12.0679645696,  23347.9184925760, -1.2754553E-04,  5.8380585E-08, -2.3407289E-10,
                7.4345550,  6.4565995078,  -1847.7052474301, -1.5181570E-04, -1.8213063E-07,  9.1183272E-10,
               -6.7314363, -4.0265854988, -16133.8556271171, -9.0955337E-05, -2.4484477E-07,  1.1116826E-09,
                6.5795750, 16.8104074692,  14323.3509980023, -2.2066770E-04, -1.1756732E-07,  5.4866364E-10,
               -6.4600721,  1.5847795630,   9061.7681128890, -6.6685176E-05, -4.3335556E-09, -3.4222998E-11,
               -6.2964773,  4.8837157343,  25300.3984729215, -1.9206388E-04, -1.4849843E-08,  6.0650192E-11,
               -5.6323538, -0.7707750092,    733.0766881638, -2.1899793E-04, -2.5474467E-07,  1.1521161E-09,
               -5.3683961,  6.8263720663,  16204.8433027325, -9.7115356E-05,  2.7023515E-08, -1.3414795E-10,
               -5.3112784,  3.9403341353,  17390.4595376141,  8.5627574E-05,  2.4607756E-07, -1.2205621E-09,
               -5.0759179,  0.6845236457,    523.5272223331,  2.1367016E-04,  2.5597745E-07, -1.2609955E-09,
               -4.8396143, -1.6710309265,  -7805.1642023920,  6.1357413E-05,  5.5663398E-09, -7.4656459E-11,
               -4.8057401,  3.5705615768,   -662.0890125485,  3.0927234E-05,  3.6923410E-08, -1.7458141E-10,
                3.9840545,  8.6945689615,  33419.5404318159,  3.9291696E-04,  7.4628340E-07, -3.5388005E-09,
                3.6744619, 19.1659620415,  22652.0424227274, -6.8354947E-05,  1.3284380E-07, -6.3767543E-10,
                2.9984815, 20.0662179587,  31190.2833132833, -3.4871029E-04, -1.2746721E-07,  5.8909710E-10,
                2.7986413, -2.5281611620, -16971.7070481963,  3.4437664E-04,  2.6526096E-07, -1.2469893E-09,
                2.4138774, 17.7106633865,  22861.5918885581, -5.0102304E-04, -3.7787833E-07,  1.7754362E-09,
                2.1863132,  5.5132179088,  -9757.6441827375,  1.2587576E-04,  7.8796768E-08, -3.6937954E-10,
                2.1461692, 13.4801375428,  23766.6709819937,  3.0245868E-04,  5.6971910E-07, -2.7016242E-09,
                1.7659832, 11.1677086523,  14809.6776020201,  1.5280981E-04,  3.1869159E-07, -1.4608454E-09,
               -1.6244212,  7.3137297434,   7318.8375983742, -4.3483492E-04, -4.4182525E-07,  2.0841655E-09,
                1.5813036,  5.4387584720,  16552.6081165349,  5.2095955E-04,  7.5618329E-07, -3.5792340E-09,
                1.5197528, 16.7359480324,  40633.6032972747,  1.7441609E-04,  5.5981921E-07, -2.6611908E-09,
                1.5156341,  1.7023646816, -17876.7861416319, -4.5910508E-04, -6.8233647E-07,  3.2300712E-09,
                1.5102092,  5.4977296450,   8399.6847301375, -3.3094061E-05,  3.1973462E-08, -1.5436468E-10,
               -1.3178223,  9.6261586339,  16275.8309783478, -2.8518605E-04, -1.9079775E-07,  8.4338673E-10,
               -1.2642739, 11.9817132061,  24604.5224030729, -1.3287330E-04,  5.9613369E-08, -3.4295235E-10,
                1.1918723, 22.4217725310,  39518.9747380084, -1.9639754E-04,  1.2294390E-07, -5.9724197E-10,
                1.1346110, 14.4235191419,  31676.6099173011,  2.4767216E-05,  3.0879170E-07, -1.4204120E-09,
                1.0857810,  8.8552797618,   5852.6842220465,  3.1609367E-06,  6.7664088E-08, -2.2006663E-10,
               -1.0193852,  7.2392703065,  33629.0898976466, -3.9751134E-05,  2.3556127E-07, -1.1256889E-09,
               -0.8227141, 11.0814572888,  16066.2815125171,  1.4748204E-04,  3.1992438E-07, -1.5697249E-09,
                0.8042238,  3.5274358950,    -33.7870573000,  2.8263353E-05,  3.7539802E-08, -2.2902113E-10,
                0.8025939,  6.7832463846,  16833.1452579809, -9.9779237E-05,  2.7639907E-08, -1.8858767E-10,
               -0.7931866, -6.3821400710, -24462.5470518423, -2.4326809E-04, -4.9525589E-07,  2.2980217E-09,
               -0.7910153,  6.3703481443,   -591.1013369332, -1.5714346E-04, -1.8089785E-07,  8.0295327E-10,
               -0.6674056,  9.1819266386,  24533.5347274576,  5.5197395E-05,  2.7743463E-07, -1.3204870E-09,
                0.6502226,  4.1010449356, -10176.3966721553, -3.0412845E-04, -4.3254175E-07,  2.0981718E-09,
               -0.6388131,  6.2958887075,  25719.1509623392,  2.3794032E-04,  4.9648867E-07, -2.4069012E-09
        };

        template <class Dummy>
        const elp2000_coefficient_t impl<Dummy>::M21[] = {
            0.0743000, 11.9537467337,  6480.9861772950,  4.9705523E-07,  6.8280480E-08, -2.7450635E-10,
            0.0304300,  8.7259027166,  7737.5900877920, -4.8307078E-06,  6.9513264E-08, -3.8338581E-10,
            0.0222900, 12.8540026510, 15019.2270678508, -2.7985829E-04, -1.9203053E-07,  9.5226618E-10,
            0.0199900, 15.2095572232, 23347.9184925760, -1.2754553E-04,  5.8380585E-08, -2.3407289E-10,
            0.0186900,  9.5981921614, -1847.7052474301, -1.5181570E-04, -1.8213063E-07,  9.1183272E-10,
            0.0169600,  7.1681781524, 16133.8556271171,  9.0955337E-05,  2.4484477E-07, -1.1116826E-09,
            0.0162300,  1.5847795630,  9061.7681128890, -6.6685176E-05, -4.3335556E-09, -3.4222998E-11,
            0.0141900, -0.7707750092,   733.0766881638, -2.1899793E-04, -2.5474467E-07,  1.1521161E-09
        };

        template <class Dummy>
        REAL impl<Dummy>::calc_moon_latitude(REAL t) {
            REAL L0 = periodic_terms(M20, sizeof(M20) / sizeof(*M20), t);
            REAL L1 = periodic_terms(M21, sizeof(M21) / sizeof(*M21), t);

            REAL L = L0 + L1 * t;
            L *= (RADIAN_PER_DEGREE / 3600);

            return L;
        }

        // 岁差
        template <class Dummy>
        REAL impl<Dummy>::adjust_precession(REAL jd, REAL l) {
            static const REAL P[] = { 50287.92262, 111.24406, 0.07699, -0.23479, -0.00178, 0.00018, 0.00001 };

            REAL t = jd / 365250.0;
            REAL t0 = 1, v = 0;
            for (auto i : P) {
                t0 *= t;
                v += i * t0;
            }

            return clamp_randians(l + (v + 2.9965 * t) * (RADIAN_PER_DEGREE / 3600));
        }

        template <class Dummy>
        REAL impl<Dummy>::calc_moon_ecliptic_longitude(REAL jd) {
            REAL t = jd / 36525.0;

            REAL l = calc_moon_longitude(t);
            l = adjust_precession(jd, l);
            return l;
        }

        template <class Dummy>
        REAL impl<Dummy>::get_moon_ecliptic_longitude(REAL jd) {
            return calc_moon_ecliptic_longitude(jd - JD2000) * DEGREE_PER_RADIAN;
        }

        // 太阳的地心黄经
        template <class Dummy>
        REAL impl<Dummy>::get_sun_ecliptic_longitude(REAL jd) {
            return calc_sun_position(jd - JD2000).longitude * DEGREE_PER_RADIAN;
        }

        // 世界时与原子时之差计算表
        template <class Dummy>
        const delta_time_t impl<Dummy>::D[] = {
            { -4000, 108371.7, -13036.80, 392.000,  0.0000 },
            {  -500,  17201.0,   -627.82,  16.170, -0.3413 },
            {  -150,  12200.6,   -346.41,   5.403, -0.1593 },
            {   150,   9113.8,   -328.13,  -1.647,  0.0377 },
            {   500,   5707.5,   -391.41,   0.915,  0.3145 },
            {   900,   2203.4,   -283.45,  13.034, -0.1778 },
            {  1300,    490.1,    -57.35,   2.085, -0.0072 },
            {  1600,    120.0,     -9.81,  -1.532,  0.1403 },
            {  1700,     10.2,     -0.91,   0.510, -0.0370 },
            {  1800,     13.4,     -0.72,   0.202, -0.0193 },
            {  1830,      7.8,     -1.81,   0.416, -0.0247 },
            {  1860,      8.3,     -0.13,  -0.406,  0.0292 },
            {  1880,     -5.4,      0.32,  -0.183,  0.0173 },
            {  1900,     -2.3,      2.06,   0.169, -0.0135 },
            {  1920,     21.2,      1.69,  -0.304,  0.0167 },
            {  1940,     24.2,      1.22,  -0.064,  0.0031 },
            {  1960,     33.2,      0.51,   0.231, -0.0109 },
            {  1980,     51.0,      1.29,  -0.026,  0.0032 },
            {  2000,     64.7,     -1.66,   5.224, -0.2905 },
            {  2150,    279.4,    732.95, 429.579,  0.0158 },
            {  6000 },
        };

        // 传入儒略日(JD2000起算),计算UTC与原子时的差(单位:日)
        template <class Dummy>
        REAL impl<Dummy>::calc_delta_t(REAL t) {
            const REAL y = t / 365.2425 + 2000;

            constexpr std::size_t length = sizeof(D) / sizeof(*D);
            const delta_time_t *p = D;
            for (std::size_t i = 0; i + 1 < length; ++i) {
                p = D + i;
                if (y < p[1].y) break;
            }

            const REAL t1 = (y - p->y) / (p[1].y - p->y - 0.0) * 10;
            const REAL d = p->a0 + (p->a1 + (p->a2 + p->a3 * t1) * t1) * t1;

            return d / 86400;
        }
    }

    typedef detail::impl<int> impl;

    static REAL make_julian_day(int year, int month, int day, int hour, int minute, REAL second) {
        if (month <= 2) {
            month += 12;
            --year;
        }
        int B = year / 400 - year / 100;
        REAL a = 365.25 * year;
        REAL b = 30.6001 * (month + 1);
        return floor(a) + floor(b) + B + 1720996.5 + day + hour / 24.0 + minute / 1440.0 + second / 86400.0;
    }

    static REAL calc_delta_t(REAL jd) {
        return impl::calc_delta_t(jd - astronomy::JD2000);
    }

    static void daytime_from_julian_day(REAL jd, daytime_t *p) {
        const REAL jdf = jd + 0.5;
        int a = (int)(jdf);
        REAL f = jdf - a;

        int cnd;
        if (a > 2299161) {
            cnd = (int)((a - 1867216.25) / 36524.25);
            a = a + 1 + cnd - (int)(cnd / 4);
        }
        a = a + 1524;

        int year = (int)((a - 122.1) / 365.25);
        cnd = a - (int)(365.25 * year);

        int month = (int)(cnd / 30.6001);
        p->day = cnd - (int)(month * 30.6001);

        year -= 4716;
        month -= 1;
        if (month > 12) month -= 12;
        if (month <= 2) year += 1;
        if (year < 1) year = year - 1;

        p->year = year;
        p->month = month;

        f *= 24.0;
        p->hour = (int)(f);
        f -= p->hour;

        f *= 60.0;
        p->minute = (int)(f);
        f -= p->minute;

        p->second = f * 60.0;
    }

    static inline REAL calc_moon_ecliptic_longitude(REAL jd) {
        return impl::calc_moon_ecliptic_longitude(jd);
    }

    static inline REAL get_moon_ecliptic_longitude(REAL jd) {
        return impl::get_moon_ecliptic_longitude(jd);
    }

    static inline REAL get_sun_ecliptic_longitude(REAL jd) {
        return impl::get_sun_ecliptic_longitude(jd);
    }
}

#endif

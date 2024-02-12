#include "astronomy.h"
#include <stdio.h>

static astronomy::REAL estimate_solar_term(int year, int angle) {
    int month = (angle + 105) / 30;
    if (month > 12) month -= 12;
    if (angle % 30 == 0) {
        return astronomy::make_julian_day(year, month, month < 8 ? 20 : 22, 0, 0, 0.0);
    }
    else {
        return astronomy::make_julian_day(year, month, month < 8 ? 4 : 7, 0, 0, 0.0);
    }
}

static astronomy::REAL calc_solar_term(int year, int idx) {
    constexpr static astronomy::REAL step = 0.000005;
    constexpr static astronomy::REAL step2 = step * 2;

    astronomy::REAL JD0, JD1, D, Dp;
    int angle = idx * 15;
    JD1 = estimate_solar_term(year, angle);
    do {
        JD0 = JD1;
        D = astronomy::get_sun_ecliptic_longitude(JD0);
        D = ((angle == 0) && (D > 345.0)) ? D - 360.0 : D;

        Dp = (astronomy::get_sun_ecliptic_longitude(JD0 + step) - astronomy::get_sun_ecliptic_longitude(JD0 - step)) / step2;
        JD1 = JD0 - (D - angle) / Dp;
    } while ((fabs(JD1 - JD0) > 1e-8));

    return JD0;
}

static astronomy::REAL clamp_degrees(astronomy::REAL d) {
    while (d < 0) d += 360;
    while (d > 360) d -= 360;
    return d;
}

static astronomy::REAL ecliptic_longitude_diff(astronomy::REAL jd) {
    return clamp_degrees(astronomy::get_moon_ecliptic_longitude(jd) - astronomy::get_sun_ecliptic_longitude(jd));
};

static astronomy::REAL estimate_new_moon_forward(astronomy::REAL jd) {
    astronomy::REAL D0, D1;
    D0 = ecliptic_longitude_diff(jd);
    for (int i = 1; i < 30; ++i) {
        jd += 1;
        D1 = ecliptic_longitude_diff(jd);
        if (D1 < D0) {
            jd -= 1;
            break;
        }
        D0 = D1;
    }
    return jd;
}

static astronomy::REAL estimate_new_moon_backward(astronomy::REAL jd) {
    constexpr astronomy::REAL ONE_DAY = 360.0 / 29.53;
    constexpr astronomy::REAL ONE_DAY_RVS = 1 / ONE_DAY;

    astronomy::REAL D0, D1;
    D0 = ecliptic_longitude_diff(jd);

    if (D0 > ONE_DAY) {
        jd -= D0 * ONE_DAY_RVS;
        D1 = ecliptic_longitude_diff(jd);
        if (D1 > D0) {
            do {
                jd += 1;
                D1 = ecliptic_longitude_diff(jd);
            } while (D1 > D0);
            return jd - 1;
        }
        else if (D1 < D0) {
            D0 = D1;
            do {
                jd -= 1;
                D1 = ecliptic_longitude_diff(jd);
            } while (D1 < D0);
            return jd;
        }

        return jd;
    }

    for (int i = 1; i < 30; ++i) {
        jd -= 1;
        D1 = ecliptic_longitude_diff(jd);
        if (D1 > D0) {
            break;
        }
        D0 = D1;
    }

    return jd;
}

static astronomy::REAL calc_new_moon_nearby(astronomy::REAL jd) {
    constexpr static astronomy::REAL step = 0.000005;
    constexpr static astronomy::REAL step2 = step * 2;

    astronomy::REAL JD0, JD1, D, Dp;
    JD1 = jd;
    do {
        JD0 = JD1;
        D = ecliptic_longitude_diff(JD0);
        D = (D > 345.0) ? D - 360.0 : D;

        Dp = (ecliptic_longitude_diff(JD0 + step) - ecliptic_longitude_diff(JD0 - step)) / step2;
        JD1 = JD0 - D / Dp;
    } while ((fabs(JD1 - JD0) > 1e-8));

    return JD0;
}

static astronomy::REAL calc_new_moon_nearby(int year, int month, int day) {
    return calc_new_moon_nearby(astronomy::make_julian_day(year, month, day, 0, 0, 0));
}

static constexpr const char *solar_terms_names[] = {
    "小寒", "大寒", "立春", "雨水", "驚蟄", "春分", "清明", "穀雨", "立夏", "小滿", "芒種", "夏至",
    "小暑", "大暑", "立秋", "處暑", "白露", "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至"
};

static constexpr const char *month_names[] = {
    "正月", "二月", "三月", "四月", "五月", "六月", "七月", "八月", "九月", "十月", "冬月", "臘月",
};

static constexpr const char *day_names[] = {
    "初一", "初二", "初三", "初四", "初五", "初六", "初七", "初八", "初九", "初十",
    "十一", "十二", "十三", "十四", "十五", "十六", "十七", "十八", "十九", "二十",
    "廿一", "廿二", "廿三", "廿四", "廿五", "廿六", "廿七", "廿八", "廿九", "三十",
};

static constexpr const char *celestial_stems[10] = {
    "甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"
};

static constexpr const char *terrestrial_branches[12] = {
    "子", "丑", "寅", "卯", "辰", "巳", "午", "未", "申", "酉", "戌", "亥"
};

// NOTE: 1928年及之前的农历用北京地方时116°23′E，1929年开始使用120°E平太阳时
// (116+23/60)*4*60=(465+8/15)*60=27932
static constexpr astronomy::REAL TIMEZONE_BEIJING = 8.0 / 24.0;
static constexpr astronomy::REAL TIMEZONE_BEIJING_LOCAL = 27932.0 / 86400.0;

static void print_daytime(const astronomy::daytime_t &dt) {
#if 1
    if ((dt.hour != 0 || dt.minute > 30) && (dt.hour != 23 || dt.minute < 30)) {
        printf("%.2d-%.2d %.2d:%.2d:%06.3f  ", dt.month, dt.day, dt.hour, dt.minute, dt.second);
    }
    else {
        printf("%.2d-%.2d %.2d:%.2d:%06.3f *", dt.month, dt.day, dt.hour, dt.minute, dt.second);
    }
#else
    printf("%.2d-%.2d", dt.month, dt.day);
#endif
}

static void print_daytime_cstb(const astronomy::daytime_t &dt) {
    int y = dt.year, m = dt.month, d = dt.day;
    if (m == 1 || m == 2) {
        m += 12;
        --y;
    }
    int n = (y * 5 + (y >> 2) - y / 100 + y / 400 + ((m - 1) & 1) * 30 + (((m - 2) * 19) >> 5) + d + 8) % 60;
    printf("%s%s", celestial_stems[n % 10], terrestrial_branches[n % 12]);
}

static int days_offset(const astronomy::daytime_t &dt) {
    int y = dt.year, m = dt.month, d = dt.day;
    if (m < 3) {
        m += 12;
        --y;
    }
    return y * 365 + (y >> 2) - y / 100 + y / 400 + (m - 3) * 30 + (((m - 2) * 19) >> 5) + d;
}

// NOTE: 一种朴素的想法，直接计算0点与24点，如果这两个时刻的值会跳转，说明节气、朔在该日
// 然而，julian_day 是有偏差的，无法反算
static void calc_solar_term_for_year(int y) {
    const astronomy::REAL tz = y >= 1929 ? TIMEZONE_BEIJING : TIMEZONE_BEIJING_LOCAL;
    astronomy::daytime_t dt;

    printf("// %.2d :", y % 100);
    for (int i = 0; i < 24; ++i) {
        astronomy::REAL jd = calc_solar_term(y, i >= 5 ? i - 5 : i + 19) + tz;
        astronomy::daytime_from_julian_day(jd - astronomy::calc_delta_t(jd), &dt);
        printf(" %d", dt.day);
    }
    printf("\n");
}

static void calc_solar_term_for_year_full(int y) {
    printf("// %.2d :\n", y % 100);
    const astronomy::REAL tz = y >= 1929 ? TIMEZONE_BEIJING : TIMEZONE_BEIJING_LOCAL;
    astronomy::daytime_t dt;

    for (int i = 0; i < 24; ++i) {
        astronomy::REAL jd = calc_solar_term(y, i >= 5 ? i - 5 : i + 19) + tz;
        astronomy::daytime_from_julian_day(jd - astronomy::calc_delta_t(jd), &dt);

        printf("%s : ", solar_terms_names[i]);
        print_daytime(dt);
        printf("\n");
    }
    printf("\n\n");
}

static void calc_new_moon_for_year_full(int y) {
    const astronomy::REAL tz = y >= 1929 ? TIMEZONE_BEIJING : TIMEZONE_BEIJING_LOCAL;
    astronomy::daytime_t dt;

    astronomy::REAL jd = astronomy::make_julian_day(y, 1, 1, 0, 0, 0.0) + tz;
    jd = calc_new_moon_nearby(estimate_new_moon_forward(jd));

    astronomy::daytime_from_julian_day(jd + tz - astronomy::calc_delta_t(jd + tz), &dt);
    printf("%.2d-%.2d %.2d:%.2d:%06.3f\n", dt.month, dt.day, dt.hour, dt.minute, dt.second);

    unsigned bit = dt.day << 12;
    int offset = days_offset(dt);

    for (int i = 0; i < 13; ++i) {
        jd = calc_new_moon_nearby(jd + 29.53);

        astronomy::daytime_from_julian_day(jd + tz - astronomy::calc_delta_t(jd + tz), &dt);
        printf("%.2d-%.2d %.2d:%.2d:%06.3f\n", dt.month, dt.day, dt.hour, dt.minute, dt.second);

        if (dt.year == y) {
            int o = days_offset(dt);
            if (o - offset == 30) {
                bit |= 1 << i;
            }
            offset = o;
        }
    }

    printf("0x%05x\n", bit);
}

#define DISPLAY_AS_CSTB 1

static void calc_chn_cal(int y) {
    constexpr int WINTER_SOLSTICE_INDEX = 23 - 5;
    const astronomy::REAL tz = y >= 1929 ? TIMEZONE_BEIJING : TIMEZONE_BEIJING_LOCAL;

    struct MyDayTime {
        astronomy::daytime_t dt;
        int ofst;

        void set(astronomy::REAL jd) {
            astronomy::daytime_from_julian_day(jd - astronomy::calc_delta_t(jd), &dt);
            ofst = days_offset(dt);
        }
    };


    // 由于农历的置闰是以冬至为锚点的，11、12月是否闰取决于上一个周期，而1~10月是否闰取决于下一个周期
    // 这里为了显示，把节气也显示出来，所以需要24*2，多出来的3是上一年的小雪、大雪、冬至
    // 朔日需要本来只需要计算26个，又因为如果冬至离朔日很近的时候，可能迭代到上一个月的，加之腊月需要显示大小，故有28
    MyDayTime solar_terms[51]{}, new_moons[28]{};

    // 上年冬至、以及上年冬至之前的朔
    astronomy::REAL jd_st = calc_solar_term(y - 1, WINTER_SOLSTICE_INDEX);
    astronomy::REAL jd_nm = calc_new_moon_nearby(estimate_new_moon_backward(jd_st));
    int st_idx = 2, nm_idx = 1;


    // 下标0和1是小雪、大雪，这两个有可能跟冬至在同一个月（概率较小）
    // 下标0预留冬至之前的朔
    solar_terms[2].set(jd_st + tz);
    new_moons[1].set(jd_nm + tz);

    // 如果朔比冬至大，则说明迭代到下一个月的朔了，需要检查更早一个朔
    if (new_moons[1].ofst > solar_terms[2].ofst) {
        astronomy::REAL jd_tmp = calc_new_moon_nearby(jd_nm - 29.53);
        new_moons[0].set(jd_tmp + tz);
        if (new_moons[0].ofst < solar_terms[2].ofst) {
            nm_idx = 0;
        }
    }
    else {
        astronomy::REAL jd_tmp = calc_new_moon_nearby(jd_nm + 29.53);
        new_moons[2].set(jd_tmp + tz);
        if (new_moons[2].ofst == solar_terms[2].ofst) {
            new_moons[1] = new_moons[2];
            jd_nm = jd_tmp;
        }
    }

    // 上年小雪、大雪
    for (int i = 0; i < 2; ++i) {
        jd_st = calc_solar_term(y - 1, (WINTER_SOLSTICE_INDEX + 22 + i) % 24);
        solar_terms[i].set(jd_st + tz);
    }

    // 上年冬至~今年冬至
    for (int i = 0; i < 24; ++i) {
        jd_st = calc_solar_term(y, i >= 5 ? i - 5 : i + 19);
        solar_terms[i + 3].set(jd_st + tz);
    }

    // 今年冬至~次年冬至
    for (int i = 0; i < 24; ++i) {
        jd_st = calc_solar_term(y + 1, i >= 5 ? i - 5 : i + 19);
        solar_terms[i + 27].set(jd_st + tz);
    }

    // 朔
    for (int i = 2; i < 28; ++i) {
        jd_nm = calc_new_moon_nearby(jd_nm + 29.53);
        new_moons[i].set(jd_nm + tz);
    }

    printf("%d\n", y);

#if 0
    printf("solar terms:\n");
    for (const auto &st : solar_terms) {
        print_daytime(st.dt);
        printf("\n");
    }

    printf("new moons:\n");
    for (int i = nm_idx; i < 28; ++i) {
        print_daytime(new_moons[i].dt);
        printf("\n");
    }
#endif

    int leap = 0;

    // 闰月在上年冬至~今年冬至区间
    // solar_terms[26]为今年冬至
    if (solar_terms[26].ofst >= new_moons[nm_idx + 13].ofst) {
        int ms_idx = 2, o;
        for (int i = 0; i < 13; ++i) {
            o = solar_terms[ms_idx].ofst;
            if (new_moons[nm_idx + i + 1].ofst <= o) {
                leap = nm_idx + i;
                break;
            }
            ms_idx += 2;
        }
    }
    // 闰月在今年冬至~下年冬至区间
    // solar_terms[50]为下年冬至
    else if (solar_terms[50].ofst >= new_moons[nm_idx + 25].ofst) {
        int ms_idx = 26, o;
        for (int i = 0; i < 13; ++i) {
            o = solar_terms[ms_idx].ofst;
            if (new_moons[nm_idx + i + 13].ofst <= o) {
                leap = nm_idx + i + 12;
                break;
            }
            ms_idx += 2;
        }
    }

    st_idx = 0;

    while (nm_idx < 16) {
        const auto &mn0 = new_moons[nm_idx];
        const auto &mn1 = new_moons[nm_idx + 1];
        bool major = (mn1.ofst - mn0.ofst) == 30;

        int readable;
        if (leap == 0 || nm_idx < leap) {
            // 无闰、或本月在闰月之前
            readable = (nm_idx + 9) % 12 + 1;
        }
        else {
            // 有闰且本月在闰月或者之后
            readable = (nm_idx + 8) % 12 + 1;
        }

#if 0
        // 按岁显示
        if (nm_idx > (leap == 0 || nm_idx <= leap ? 13 : 14)) {
            break;
        }
#endif

#if 1
        // 按年显示
        if (leap == 0 ? (nm_idx < 3 || nm_idx > 14) : (nm_idx < (leap >= 4 ? 3 : 4) || nm_idx > (leap <= 15 ? 15 : 14))) {
            ++nm_idx;
            continue;
        }
#endif


#if DISPLAY_AS_CSTB
        // 汉字方式显示月份 天支日期
        printf("%s%s%s ", leap == 0 || nm_idx != leap ? "　" : "閏", month_names[readable - 1], major ? "大" : "小");
        print_daytime_cstb(mn0.dt);
#else
#if 0
        // 数字方式显示月份
        printf("%c", leap == 0 || nm_idx != leap ? ' ' : '+');
        printf("%.2d ", readable);
        printf("%c (", major ? '+' : '-');
        print_daytime(mn0.dt);
        printf(")");
#endif

#if 0
        // 汉字方式显示月份
        printf("%s%s%s (", leap == 0 || nm_idx != leap ? "　" : "閏", month_names[readable - 1], major ? "大" : "小");
        print_daytime(mn0.dt);
        printf(")");
#endif
#endif

        ++nm_idx;

        // 显示节气
        while (1) {
            const auto &st = solar_terms[st_idx];

            // 节气超过本月了
            if (st.ofst >= mn1.ofst) {
                break;
            }

            // 节气落在本月内
            if (st.ofst >= mn0.ofst) {
#if DISPLAY_AS_CSTB
                printf(" %s", day_names[st.ofst - mn0.ofst]);
                print_daytime_cstb(st.dt);
                printf("%s", solar_terms_names[(st_idx + 21) % 24]);
#else
                printf(" %s (", solar_terms_names[(st_idx + 21) % 24]);
                print_daytime(st.dt);
                printf(")");
#endif
            }

            ++st_idx;
        }

        printf("\n");
    }

}

int main() {
    //calc_new_moon_for_year_full(2024);
    //calc_solar_term_for_year_full(2022);

    //calc_solar_term_for_year_full(1900);
    //calc_solar_term_for_year_full(1912);
    //calc_solar_term_for_year_full(1913);
    //calc_solar_term_for_year_full(1917);
    //calc_solar_term_for_year_full(1923);
    //calc_solar_term_for_year_full(1927);
    //calc_solar_term_for_year_full(1928);
    //calc_solar_term_for_year_full(1844);
    //calc_solar_term_for_year_full(1846);
    //calc_solar_term_for_year_full(1849);
    //calc_solar_term_for_year_full(1850);
    //calc_solar_term_for_year_full(1851);
    //calc_solar_term_for_year_full(1855);
    //calc_solar_term_for_year_full(1862);
    //calc_solar_term_for_year_full(1864);
    //calc_solar_term_for_year_full(1866);
    //calc_solar_term_for_year_full(1867);
    //calc_solar_term_for_year_full(1879);
    //calc_solar_term_for_year_full(1883);
    //calc_solar_term_for_year_full(1884);
    //calc_solar_term_for_year_full(1886);
    //calc_solar_term_for_year_full(1895);
    //calc_solar_term_for_year_full(1898);
    //calc_solar_term_for_year_full(1899);

    //printf("===\n");
    //for (int i = 1900; i < 2000; ++i) {
    //    calc_solar_term_for_year(i);
    //}
    //for (int i = 1899; i >= 1799; --i) {
    //    calc_solar_term_for_year(i);
    //}
    //calc_solar_term_for_year_full(2023);
    //calc_solar_term_for_year_full(2024);
    //calc_solar_term_for_year_full(2025);

    //calc_new_moon_for_year_full(2032);
    //calc_chn_cal(2033);
    //calc_chn_cal(2034);

#if 1
    // 测试数据2262年 闰正月
    calc_chn_cal(2261);
    calc_chn_cal(2262);
    calc_chn_cal(2263);

    // 2023 闰二月
    calc_chn_cal(2022);
    calc_chn_cal(2023);
    calc_chn_cal(2024);

    // 1993 闰三月
    calc_chn_cal(1992);
    calc_chn_cal(1993);
    calc_chn_cal(1994);

    // 2020 闰四月
    calc_chn_cal(2019);
    calc_chn_cal(2020);
    calc_chn_cal(2021);

    // 2009 闰五月
    calc_chn_cal(2008);
    calc_chn_cal(2009);
    calc_chn_cal(2010);

    // 2017 闰六月
    calc_chn_cal(2016);
    calc_chn_cal(2017);
    calc_chn_cal(2018);

    // 2006 闰七月
    calc_chn_cal(2005);
    calc_chn_cal(2006);
    calc_chn_cal(2007);

    // 1995 闰八月
    calc_chn_cal(1994);
    calc_chn_cal(1995);
    calc_chn_cal(1996);

    // 2014 闰九月
    calc_chn_cal(2013);
    calc_chn_cal(2014);
    calc_chn_cal(2015);

    // 1984 闰十月
    calc_chn_cal(1983);
    calc_chn_cal(1984);
    calc_chn_cal(1985);

    // 2033 闰十一月
    calc_chn_cal(2032);
    calc_chn_cal(2033);
    calc_chn_cal(2034);

    // 测试数据3358年 闰十二月
    calc_chn_cal(3357);
    calc_chn_cal(3358);
    calc_chn_cal(3359);
#endif

#if 0
    //         计算                历书
    // 1803-09-08 23:52:11.667 | 09 白露
    // 1805-08-23 23:51:55.314 | 24 处暑 七月大
    // 1807-02-05 00:00:17.045 | 04 立春
    // 1808-06-21 23:51:45.160 | 22 夏至
    // 1809-11-22 23:44:42.817 | 23 小雪
    // 1812-03-06 00:19:51.796 | 05 惊蛰
    // 1813-04-30 23:56:19.686 | 01 四月小
    // 1814-12-07 23:58:25.509 | 08 大雪
    // 1817-10-08 23:40:19.648 | 09 寒露
    //      10-10 23:47:39.041 | 11 九月小
    // 1818-09-23 23:54:52.482 | 24 秋分
    // 1820-02-20 00:06:22.697 | 19 雨水
    //      12-05 23:47:16.656 | 06 十一月小
    // 1823-05-10 23:55:25.652 | 11 四月小 
    // 1824-08-07 23:44:18.970 | 08 立秋
    // 1825-11-07 23:51:41.744 | 08 立冬
    // 1829-10-23 23:55:06.218 | 24 霜降
    //      11-07 23:21:08.805 | 08 立冬
    // 1836-09-07 23:37:15.542 | 08 白露
    // 1839-01-20 23:59:09.959 | 21 大寒
    // 1842-01-05 23:55:49.309 | 06 小寒
    //      01-12 00:01:08.204 | 11 十二月大
    //      11-02 23:54:28.000 | 03 十月小
    //      11-22 23:56:57.039 | 23 小雪
    // 1844-06-05 23:35:02.385 | 06 芒种
    // 1846-11-22 23:09:36.239 | 23 小雪
    // 1848-12-21 23:45:21.945 | 22 冬至
    // 1849-09-16 23:47:33.521 | 17 八月小
    // 1850-10-08 23:23:59.024 | 09 寒露
    // 1851-09-23 23:36:00.653 | 24 秋分
    //      12-07 23:29:24.869 | 08 大雪
    // 1856-11-27 23:46:58.303 | 28 十一月小
    // 1861-11-02 23:50:25.224 | 03 十月小
    // 1862-10-23 23:30:45.572 | 24 霜降
    //      11-07 23:09:11.430 | 08 立冬
    // 1864-07-22 23:36:01.878 | 23 大暑
    // 1865-09-07 23:52:34.285 | 08 白露
    // 1866-10-23 22:59:36.710 | 24 霜降
    // 1867-07-07 23:36:24.529 | 08 小暑
    //      08-23 23:38:56.617 | 24 处暑
    // 1869-05-11 23:52:26.505 | 12 四月小
    // 1878-05-05 23:59:55.681 | 06 立夏
    // 1879-01-05 23:36:40.508 | 06 小寒
    //      10-08 23:50:13.701 | 09 寒露
    //      11-22 23:15:47.989 | 23 小雪
    // 1880-09-22 23:52:20.616 | 23 秋分
    //      11-02 23:41:05.343 | 03 十月小
    // 1881-12-21 23:46:13.515 | 22 冬至
    // 1883-10-08 23:03:51.174 | 09 寒露
    // 1884-09-22 23:06:57.334 | 23 秋分
    //      12-06 23:35:51.495 | 07 大雪
    // 1886-08-07 23:29:50.334 | 08 立秋
    // 1887-03-24 23:54:58.213 | 25 三月小
    // 1893-07-22 23:51:14.096 | 23 大暑
    // 1895-10-23 23:32:12.762 | 24 霜降
    //      11-07 23:22:05.606 | 08 立冬
    // 1896-07-06 23:52:05.616 | 07 小暑
    //      08-22 23:50:36.246 | 23 处暑
    // 1898-09-07 23:24:58.401 | 08 白露
    // 1899-06-21 23:31:13.699 | 22 夏至
    //      10-23 22:52:11.480 | 24 霜降
    // 1906-04-23 23:51:31.408 | 24 四月小
    // 1909-01-20 23:56:35.795 | 21 大寒
    // 1911-05-06 23:45:58.481 | 07 立夏
    // 1912-01-06 23:53:06.441 | 07 小寒
    //      10-08 23:52:18.898 | 09 寒露
    //      11-22 23:33:43.341 | 23 小雪
    // 1913-09-23 23:38:19.512 | 24 秋分
    // 1979-01-21 00:00:02.815 | 20 大寒
    // 2008-05-21 00:00:54.920 | 20 小满

    calc_chn_cal(1803);  // 09-08 23:52:11.667 | 09 白露
    calc_chn_cal(1804);  // 七月大 08-05
    calc_chn_cal(1805);  // 08-23 23:51:55.314 | 24 处暑 七月大
    calc_chn_cal(1807);  // 02-05 00:00:17.045 | 04 立春
    calc_chn_cal(1808);  // 06-21 23:51:45.160 | 22 夏至
    calc_chn_cal(1809);  // 11-22 23:44:42.817 | 23 小雪
    calc_chn_cal(1812);  // 03-06 00:19:51.796 | 05 惊蛰
    calc_chn_cal(1813);  // 04-30 23:56:19.686 | 01 四月小
    calc_chn_cal(1814);  // 12-07 23:58:25.509 | 08 大雪
    calc_chn_cal(1815);  // 清明 04-05
    calc_chn_cal(1817);  // 10-08 23:40:19.648 | 09 寒露
                         // 10-10 23:47:39.041 | 11 九月小
    calc_chn_cal(1818);  // 09-23 23:54:52.482 | 24 秋分
    calc_chn_cal(1820);  // 02-20 00:06:22.697 | 19 雨水
                         // 12-05 23:47:16.656 | 06 十一月小
    calc_chn_cal(1823);  // 05-10 23:55:25.652 | 11 四月小
    calc_chn_cal(1824);  // 08-07 23:44:18.970 | 08 立秋
    calc_chn_cal(1825);  // 11-07 23:51:41.744 | 08 立冬
    calc_chn_cal(1826);  // 小满 05-21
    calc_chn_cal(1829);  // 10-23 23:55:06.218 | 24 霜降
                         // 11-07 23:21:08.805 | 08 立冬
    calc_chn_cal(1831);  // 三月大 04-12
    calc_chn_cal(1836);  // 09-07 23:37:15.542 | 08 白露
    calc_chn_cal(1839);  // 01-20 23:59:09.959 | 21 大寒
    calc_chn_cal(1842);  // 01-12 00:01:08.204 | 11 十二月大
                         // 11-02 23:54:28.000 | 03 十月小
                         // 11-22 23:56:57.039 | 23 小雪
    calc_chn_cal(1844);  // 06-05 23:35:02.385 | 06 芒种
    calc_chn_cal(1846);  // 11-22 23:09:36.239 | 23 小雪
    calc_chn_cal(1848);  // 12-21 23:45:21.945 | 22 冬至
    calc_chn_cal(1849);  // 09-16 23:47:33.521 | 17 八月小
    calc_chn_cal(1850);  // 10-08 23:23:59.024 | 09 寒露
    calc_chn_cal(1851);  // 09-23 23:36:00.653 | 24 秋分
                         // 12-07 23:29:24.869 | 08 大雪
    calc_chn_cal(1855);  // 谷雨 04-20
    calc_chn_cal(1856);  // 11-27 23:46:58.303 | 28 十一月小
    calc_chn_cal(1861);  // 11-02 23:50:25.224 | 03 十月小
    calc_chn_cal(1862);  // 10-23 23:30:45.572 | 24 霜降
                         // 11-07 23:09:11.430 | 08 立冬
    calc_chn_cal(1863);  // 腊月大 01-19
    calc_chn_cal(1864);  // 07-22 23:36:01.878 | 23 大暑
    calc_chn_cal(1865);  // 09-07 23:52:34.285 | 08 白露
    calc_chn_cal(1866);  // 10-23 22:59:36.710 | 24 霜降
    calc_chn_cal(1867);  // 07-07 23:36:24.529 | 08 小暑
                         // 08-23 23:38:56.617 | 24 处暑
    calc_chn_cal(1869);  // 05-11 23:52:26.505 | 12 四月小
    calc_chn_cal(1878);  // 05-05 23:59:55.681 | 06 立夏
    calc_chn_cal(1879);  // 01-05 23:36:40.508 | 06 小寒
                         // 10-08 23:50:13.701 | 09 寒露
                         // 11-22 23:15:47.989 | 23 小雪
    calc_chn_cal(1880);  // 09-22 23:52:20.616 | 23 秋分
                         // 11-02 23:41:05.343 | 03 十月小
    calc_chn_cal(1881);  // 12-21 23:46:13.515 | 22 冬至
    calc_chn_cal(1883);  // 10-08 23:03:51.174 | 09 寒露
    calc_chn_cal(1884);  // 09-22 23:06:57.334 | 23 秋分
                         // 12-06 23:35:51.495 | 07 大雪
    calc_chn_cal(1886);  // 08-07 23:29:50.334 | 08 立秋
    calc_chn_cal(1887);  // 03-24 23:54:58.213 | 25 三月小
    calc_chn_cal(1893);  // 07-22 23:51:14.096 | 23 大暑
    calc_chn_cal(1895);  // 10-23 23:32:12.762 | 24 霜降
                         // 11-07 23:22:05.606 | 08 立冬
    calc_chn_cal(1896);  // 正月大 02-13
                         // 07-06 23:52:05.616 | 07 小暑
                         // 08-22 23:50:36.246 | 23 处暑
    calc_chn_cal(1898);  // 09-07 23:24:58.401 | 08 白露
    calc_chn_cal(1899);  // 06-21 23:31:13.699 | 22 夏至
                         // 10-23 22:52:11.480 | 24 霜降
    calc_chn_cal(1906);  // 04-23 23:51:31.408 | 24 四月小
    calc_chn_cal(1909);  // 01-20 23:56:35.795 | 21 大寒
    calc_chn_cal(1911);  // 05-06 23:45:58.481 | 07 立夏
    calc_chn_cal(1912);  // 01-06 23:53:06.441 | 07 小寒
                         // 10-08 23:52:18.898 | 09 寒露
                         // 11-22 23:33:43.341 | 23 小雪
    calc_chn_cal(1913);  // 09-23 23:38:19.512 | 24 秋分
    calc_chn_cal(1914);  // 十月大 11-17
    calc_chn_cal(1916);  // 正月大 02-03
    calc_chn_cal(1917);  // 大雪 12-07
    calc_chn_cal(1920);  // 十月大 11-10
    calc_chn_cal(1927);  // 白露 09-08
    calc_chn_cal(1928);  // 夏至 06-21
    calc_chn_cal(1979);  // 01-21 00:00:02.815 | 20 大寒
    calc_chn_cal(2008);  // 05-21 00:00:54.920 | 20 小满
#endif

    return 0;
}

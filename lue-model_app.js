function calculate_par(par_param, start, end, interval) {
    if (par_param == 'MODIS18C2') {
        var par = ee.ImageCollection('MODIS/061/MCD18C2')
            .filterDate(start, end)
            .map(allign_image)
            .sum()
            .reduce(ee.Reducer.sum());
        //The PAR data is updated every 3 hours, and the total daily PAR is calculated as the sum of all PAR values for that day (i.e., 0000_PAR + 0030_PAR + ... + 2100_PAR). So, after summing up all par band's value of given period, we need to multiply it with 3 * 60 * 60. (Because MODIS/061/MCD18C2 PAR band's unit is J s-1 m-2 (aka W m-2), we need to convert it to MJ 3h-1 m−2)
        var par = par.multiply(3 * 60 * 60 * 1e-6); // unit: MJ m-2 3h-1
    } else if (par_param == 'BESS') {
        var par = ee.ImageCollection('SNU/ESL/BESS/Rad/v1')
            .select(['PAR_Daily'])
            .filterDate(start, start.advance(interval + 1, 'day'))  // This dataset is not continuous, with breaks on certain dates. The interval should be more than 2 days.
            .map(allign_image)
            .sum();
        var par = par.multiply(24 * 60 * 60 * 1e-6); // Convert the unit from J s-1 m-2 to MJ d-1 m-2
    } else if (par_param == "ERA5-Land") {
        var par = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
            .select(['surface_solar_radiation_downwards_sum'])
            .filterDate(start, end)
            .map(allign_image)
            .sum()
            .multiply(0.48);  // photosynthetically active radiation (PAR) is 50% of solar shortwave radiation downwards
        var par = par.multiply(1e-6); // This is an daily-aggregated data collection. We need to convert the unit from J d-1 m-2 to MJ d-1 m-2
    } else if (par_param == "MERRA-2") {
        var par = ee.ImageCollection("NASA/GSFC/MERRA/rad/2")
            .select(['SWGNT'])
            .filterDate(start, end)
            .map(allign_image)
            .sum()
            .multiply(0.48);  // photosynthetically active radiation (PAR) is 50% of solar shortwave radiation downwards
        var par = par.multiply(60 * 60 * 1e-6); // This is an hourly time-averaged data collection. We need to convert the unit from J s-1 m-2 to MJ h-1 m-2
    } else if (par_param == "GLDAS-2.2") {
        var par = ee.ImageCollection("NASA/GLDAS/V022/CLSM/G025/DA1D")
            .select(['Swnet_tavg'])
            .filterDate(start, end)
            .map(allign_image)
            .sum()
            .multiply(0.48);  // photosynthetically active radiation (PAR) is 50% of solar shortwave radiation downwards
        var par = par.multiply(24 * 60 * 60 * 1e-6); // This is an daily-averaged (not 3-hourly-averaged) data collection. We need to convert the unit from J s-1 m-2 to MJ d-1 m-2
    } else if (par_param == "TerraClimate") {
        var par = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
            .select(['srad'])
            .filter(ee.Filter.date(start.update(null, null, 1), start.update(null, null, 2)))
            .map(allign_image)
            .sum()
            .multiply(0.1)  // scale factor is 0.1
            .multiply(0.48);  // photosynthetically active radiation (PAR) is 48% of solar shortwave radiation downwards
        var par = par.multiply(30.42 * 24 * 60 * 60 * 1e-6); // This is an monthly-averaged data collection. We need to convert the unit from J s-1 m-2 to MJ month-1 m-2
        var par = par.multiply(interval / 30.42)
    };
    return par;
}

function calculate_fpar(fpar_param, start, end, sr_max, sr_min, interval) {
    if (interval <= 8 || end.difference(start, 'day').lte(8)) {
        var end8 = start.advance(8, 'day');
    }
    // MODIS/MOD09GA_006_NDVI is not continuous, with breaks on certain dates. The interval should be at least 3 days.
    if (fpar_param == 'MOD15A2H') {
        var fpar = ee.ImageCollection('MODIS/061/MOD15A2H')
            .filterDate(start, end8)  // 2000-02-18T00:00:00Z–2023-04-23T00:00:00
            .select('Fpar_500m')
            .map(allign_image)
            .mean()
            .multiply(0.01);  // scale factor is 0.01
    } else if (fpar_param == 'NDPI') {
        var img = ee.ImageCollection("MODIS/061/MOD09A1")
            .filterDate(start, end8)  // 2000-02-18T00:00:00Z–2023-04-23T00:00:00
            .map(allign_image)
            .mean()
            .multiply(0.0001);  // scale factor is 0.0001
        var ndpi = img.expression(
            '(nir - (0.74 * red + 0.26 * swir1)) / (nir + (0.74 * red + 0.26 * swir1))',
            {
                'red': img.select('sur_refl_b01'),
                'nir': img.select('sur_refl_b02'),
                'swir1': img.select('sur_refl_b06')
            }).float()
        var fpar = ndpi  // Replace fpar with NDPI
    } else if (fpar_param == "LAI-based") {
        // fPAR = 1 - exp(-k * LAI), k is extinction coefficient (0.5, constant). Ref: Evaluation of MODIS NPP and GPP products across multiple biomes
        var lai = ee.ImageCollection("MODIS/061/MOD15A2H")
            .filterDate(start, end8)  // 2000-02-18T00:00:00Z–2023-04-23T00:00:00
            .select('Lai_500m')
            .map(allign_image)
            .mean()
            .multiply(0.1);  // scale factor is 0.1
        var fpar = lai.expression('1 - exp(-0.5 * lai)', { 'lai': lai });
    } else if (fpar_param == 'CASA') {
        // fPAR = min(SR / (SRmax - SRmin) - SRmin / (SRmax - SRmin), 0.95)
        var sr = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .map(function (image) {
                return image.expression('(1 + ndvi) / (1 - ndvi)', { 'ndvi': image });
            })
            .mean();
        var fpar = sr.expression('(sr - sr_min) / (sr_max - sr_min)', { 'sr': sr, 'sr_max': sr_max, 'sr_min': sr_min }).min(0.95);
    } else if (fpar_param == 'GLO-PEM') {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .mean()
            .multiply(1.67).subtract(0.08);
    } else if (fpar_param == 'C-Fix') {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .mean()
            .multiply(0.8642).subtract(0.0814);
    } else if (fpar_param == 'VPM') {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_EVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('EVI')
            .map(allign_image)
            .mean()
            .subtract(0.1).multiply(1.25);
    } else if (fpar_param == 'EC-LUE-V3') {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .mean()
            .multiply(1.24).subtract(0.168);
    } else if (fpar_param == 'TPM') {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .mean()
            .subtract(0.1).divide(0.8)
            .clamp(0, 1);
    } else if (fpar_param == "Wang-model") {
        var fpar = ee.ImageCollection('MODIS/MOD09GA_006_NDVI')
            .filterDate(start, end8)  // 2000-02-24T00:00:00Z–2023-02-17T00:00:00
            .select('NDVI')
            .map(allign_image)
            .mean()
            .subtract(0.05);
    };
    return fpar;
}

function calculate_lue_ci(lue_ci_param, c3c4_map, lulc, ci, ci_max, ci_min, resolution) {
    if (lue_ci_param == 'MOD17') {
        var bplut = calculate_bplut_modis(lulc);
        var lue_ci = bplut[2];
    } else if (lue_ci_param == 'P-model') {
        var Mc = 12.0107; // molecular mass of C
        var kphio = 0.081785; // Apparent quantum yield efficiency (unitless)
        var lue_ci = ee.Image(Mc * kphio);
    } else if (lue_ci_param == 'C-Fix') {
        var lue_ci = ee.Image(1.1);  // 1.1 gC/MJ (APAR)
    } else if (lue_ci_param == 'VPM') {
        // C3 = 0.42, C4 = 0.63 gC/mol APAR.
        var lue_ci = ee.Image(0);
        var value_dict = [
            { "c3c4": 10, "lue_max": 0.42 },
            { "c3c4": 20, "lue_max": 0.63 },
            { "c3c4": 30, "lue_max": 0.525 },
            { "c3c4": 40, "lue_max": 0 },
        ];
        for (var i = 0; i < value_dict.length; i++) {
            var value = value_dict[i];
            var lue_ci = lue_ci.where(c3c4_map.eq(value["c3c4"]), value["lue_max"]);
        }
    } else if (lue_ci_param == 'TEC') {
        // C3 = 1.8, C4 = 2.76 gC/mol APAR.
        var lue_ci = ee.Image(0);
        var value_dict = [
            { "c3c4": 10, "lue_max": 1.8 },
            { "c3c4": 20, "lue_max": 2.76 },
            { "c3c4": 30, "lue_max": 2.28 },
            { "c3c4": 40, "lue_max": 0 },
        ];
        for (var i = 0; i < value_dict.length; i++) {
            var value = value_dict[i];
            var lue_ci = lue_ci.where(c3c4_map.eq(value["c3c4"]), value["lue_max"]);
        }
    } else if (lue_ci_param == 'CFLUX | CI-EF') {
        /*
        f(LUE_CI) = (lue_max - lue_cs) × CInor + lue_cs;
        CInor = (CI - CImin) / (CImax - CImin)
      
        f(CI): Cloudiness index constraint (dimensionless);
        CInor: Normalized CI (unitless);
        CI: Cloudiness index (dimensionless)
        mu: mu indicates an overall sensitivity of GPP to CI, euqal to 0.46 (unitless);
        */

        var ci_nor = ci.expression('(CI - CImin) / (CImax - CImin)', { 'CI': ci, 'CImin': ci_min, 'CImax': ci_max });
        var bplut = calculate_bplut_CI_EF(lulc);
        var lue_max = bplut[2];
        var lue_cs = bplut[3];
        var lue_ci = ci_nor.expression('(lue_max - lue_cs) * CInor + lue_cs', { 'lue_max': lue_max, 'lue_cs': lue_cs, 'CInor': ci_nor });
    } else if (lue_ci_param == 'Wang-model') {
        /*
        Ref: Incorporating diffuse radiation into a light use efficiency and evapotranspiration model: An 11-year study in a high latitude deciduous forest
      
        f(LUE_CI) = lue_max * f(CI)
        f(CI) = 1 - mu * (1 - CInor)
        CInor = (CI - CImin) / (CImax - CImin)
      
        lue_max: Maximum LUE, 2.97/4.29 gC m-2 MJ-1 (without CI / with CI);
        f(CI): Cloudiness index constraint (dimensionless);
        CInor: Normalized CI (unitless);
        CI: Cloudiness index (dimensionless);
        mu: pronunce of μ, μ indicates an overall sensitivity of GPP to CI, euqal to 0.46 (unitless);
        */

        var ci_nor = ci.expression('(CI - CImin) / (CImax - CImin)', { 'CI': ci, 'CImin': ci_min, 'CImax': ci_max });
        var fci = ci_nor.expression('1 - mu * (1 - CInor)', { 'mu': 0.46, 'CInor': ci_nor });
        var lue_ci = fci.multiply(4.29);
    } else if (lue_ci_param == 'MuSyQ') {
        var bplut = calculate_bplut_MuSyQ(lulc);
        var lue_su = bplut[0];
        var lue_sh = bplut[1];
        var lue_ci = lue_su.expression('lue_su * (1 - ci) + lue_sh * ci', { 'lue_su': lue_su, 'lue_sh': lue_sh, 'ci': ci });
    } else if (lue_ci_param == 'CCW') {
        // f(CI) = 1 - K1 * CI, CI is the ratio of actual radiation to clear-sky radiation
        var bplut = calculate_bplut_ccw(lulc)
        var fci = ci.expression('1 - K1 * CI', { 'K1': bplut[2], 'CI': ci.multiply(-1).add(1) });
        var lue_ci = fci.multiply(bplut[0]);
    } else if (lue_ci_param == 'None') {
        var lue_ci = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
    };
    return lue_ci;
}

function calculate_ft(ft_param, T, Tk, lulc, Topt_CASA, c3c4_map, resolution) {
    if (ft_param == 'P-model') {
        var ft = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
        var ft_c3 = T.expression('0.352 + 0.022 * T - 0.00034 * T * T', { 'T': T });
        var ft_c4 = T.expression('-0.064 + 0.03 * T - 0.000464 * T * T', { 'T': T });
        var ft_c3c4 = ft_c3.add(ft_c4).divide(2);
        var ft = ft.where(c3c4_map.eq(10), ft_c3).where(c3c4_map.eq(20), ft_c4).where(c3c4_map.eq(30), ft_c3c4);
    } else if (ft_param == 'CFLUX | CI-EF | MOD17 | TL-LUE') {
        var bplut = calculate_bplut_modis(lulc);
        var Tmin_min = bplut[0];
        var Tmin_max = bplut[1];
        var ft = T.expression('(T - Tmin_min) / (Tmin_max - Tmin_min)', { 'T': T, 'Tmin_min': Tmin_min, 'Tmin_max': Tmin_max }).clamp(0, 1);
    } else if (ft_param == 'CASA') {
        /*
        f(T) = Tε1 * Tε2
        Tε1 = 0.8 + 0.02 * Topt - 0.0005* Topt ** 2
        Tε2 = 1.1814 / (1 + exp(0.2 * (Topt - 10 -T))) / (1 + exp(0.3 * (-Topt - 10 + T)))
        
        In this equation, Topt represents the optimal temperature for plant growth, which is defined as the monthly average temperature (℃) when the NDVI value reaches its maximum in a certain region within a year. When the average temperature in a month is less than or equal to -10℃, Tε1 is set to 0.
        When the monthly average temperature T is 10℃ higher or 13℃ lower than the optimal temperature Topt, the value of Tε2 for that month is equal to half of the value of Tε2 when the monthly average temperature T is equal to the optimal temperature Topt, which is 0.4956 [5].
    
        Due to a typo in the original Potter formula [1], several papers were misled [2-3]. Here, we quote the formula in Field CB 1995's paper[4].
        Ref: 
          1. Potter, Christopher S., et al. "Terrestrial ecosystem production: a process model based on global satellite and surface data." Global biogeochemical cycles 7.4 (1993): 811-841.
          2. Jiang, Shouzheng, et al. "Comparison of satellite-based models for estimating gross primary productivity in agroecosystems." Agricultural and Forest Meteorology 297 (2021): 108253.
          3. Bao, Shanning, et al. "Environment-sensitivity functions for gross primary productivity in light use efficiency models." Agricultural and Forest Meteorology 312 (2022): 108708.
          4. Field, Christopher B., James T. Randerson, and Carolyn M. Malmström. "Global net primary production: combining ecology and remote sensing." Remote sensing of Environment 51.1 (1995): 74-88.
          5. 朱文泉, 潘耀忠, and 张锦水. "中国陆地植被净初级生产力遥感估算." 植物生态学报 31.3 (2007): 413-424.
        */

        var T1 = T.expression('0.8 + 0.02 * Topt - 0.0005 * Topt * Topt', { 'Topt': Topt_CASA });
        var T2 = T.expression('1.1814 / (1 + exp(0.2 * (Topt - 10 - T))) / (1 + exp(0.3 * (-Topt - 10 + T)))', { 'T': T, 'Topt': Topt_CASA });
        var T1 = T1.where(T1.lte(-10), 0);
        var T2 = T2.where(T.gt(Topt_CASA.add(10)).or(T.lt(Topt_CASA.subtract(13))), 0.4956);
        var ft = T1.multiply(T2);
    } else if (ft_param == 'C-Fix') {
        /*
        Ref: Estimation of carbon mass fluxes over Europe using the C-Fix model and Euroflux data, section 2.1.2.
    
        f(T) = exp(C1 - ΔHa / Rg / T) / (1 + exp((ΔS * T - ΔHd) / Rg / T))
        
        f(T): The temperature dependency factor.
        C1: A constant, equal to 21.77.
        ΔHa: The activation energy of photosynthesis, equal to 52750 J mol-1.
        T: The air temperature in Kelvin.
        ΔS: The entropy of the denaturationequilibrium of CO, equal to 704.95 J K-1 mol -1.
        ΔHd: The deactivation energy of photosynthesis, equal to 211000 J mol-1.
        Rg: The gas constant, equal to 8.314 J mol-1 K-1.
        */

        var ft = Tk.expression('exp(C1 - DeltaHa / Rg / T) / (1 + exp((DeltaS * T - DeltaHd) / Rg / T))', { 'C1': 21.77, 'DeltaHa': 52750, 'T': Tk, 'DeltaS': 704.95, 'DeltaHd': 211000, 'Rg': 8.314 });
    } else if (ft_param == 'CCW | EC-LUE-V3 | GLO-PEM | MVPM | TEC | VPM') {
        var expr = '((T - Tmin) * (T - Tmax)) / ((T - Tmin) * (T - Tmax) - (T - Topt) * (T - Topt))'
        var ft = T.expression('((T - Tmin) * (T - Tmax)) / ((T - Tmin) * (T - Tmax) - (T - Topt) * (T - Topt))', { 'T': T, 'Tmax': 40, 'Tmin': 0, 'Topt': 20.33 }).clamp(0, 1);
    } else if (ft_param == 'CI-LUE') {
        var ft = T.expression('(T - Tmin) / (Tmax - Tmin)', { 'T': T, 'Tmax': 40, 'Tmin': 0 }).clamp(0, 1);
    } else if (ft_param == 'MuSyQ') {
        /*
        Ref: New Global MuSyQ GPP/NPP Remote Sensing Products From 1981 to 2018.
        Original paper 'Estimating Vegetation Primary Production in the Heihe River Basin of China with Multi-Source and Multi-Scale Data' (Cui et al, 2016b) use different version of formular, which is slightly more complicate.
        */

        var expr = '1 / ((1 + exp(0.2 * (Topt - 10 - T))) * (1 + exp(0.3 * (-Topt - 10 + T))))'
        var Topt = calculate_bplut_MuSyQ(lulc)[2];
        var ft = T.expression(expr, { 'T': T, 'Topt': Topt });
    } else if (ft_param == 'Wang-model') {
        /*
        Ref: Incorporating diffuse radiation into a light use efficiency and evapotranspiration model: An 11-year study in a high latitude deciduous forest.
    
        f(T) = 1.1814 / (1 + exp(-To - 10 + T)) / (1 + exp(To - 10 - T))
    
        f(T): The temperature dependency factor.
        To: The optimum plant growth temperature, Ta at max{PAR·fAPAR·Ta/VPD}, 16.51 °C for this study
        T: The air temperature in °C.
        */

        var ft = T.expression('1.1814 / (1 + exp(-To - 10 + T)) / (1 + exp(To - 10 - T))', { 'T': T, 'To': 16.51 });
    } else if (ft_param == 'None') {
        var ft = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
    };
    return ft;
}

function calculate_fw(fw_param, start, end, T, Tdew, Pmst, sm, sm_max, sm_min, lswi_max, actual_eva, potential_eva, resolution) {
    if (fw_param == 'P-model') {
        /*
        Ref: https://pyrealm.readthedocs.io/en/latest/users/pmodel/pmodel_details/soil_moisture.html; https://cran.r-project.org/web/packages/rpmodel/rpmodel.pdf
        
        β = q * (θ − θ∗) ** 2 + 1
        q = (a + bα − 1) / (θ∗ − θ0) ** 2
    
        β: Soil moisture stress factor;
        θ: Relative soil moisture as a fraction of field capacity;
        θ∗: An upper bound in relative soil moisture, set to 1.0;
        θ0: A lower bound in relative soil moisture, set to 0.0;
        q: aridity sensitivity parameter;
        a, b: Parameter  determining  the sensitivity of the empirical soil moisture stress function. Defaults to 0.0, 0.685;
        α: Local monthly mean (actual evapotranspiration / potential evapotranspiration), measure for average aridity;
        
        --> β = 1 - (1 - 0.685α) * (1 - θ) ** 2 <--
        */

        var start = start.update(null, null, 1);
        var end = start.advance(1, 'month');
        var evaporation_month_mean = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
            .filter(ee.Filter.date(start, end))
            .select(['evaporation_from_vegetation_transpiration_sum', 'potential_evaporation_sum'])
            .map(allign_image)
            .mean();
        var alpha = evaporation_month_mean.select('evaporation_from_vegetation_transpiration_sum').divide(evaporation_month_mean.select('potential_evaporation_sum'));
        var fw = alpha.expression('1 - (1 - 0.685 * alpha) * (1 - sm) ** 2', { 'alpha': alpha, 'sm': sm });
    } else if (fw_param == 'GLO-PEM') {
        /*
        Ref: https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html; https://www.sciencedirect.com/science/article/pii/S0304380015000137
    
        f(W) = f(SHD) * f(SM)
        f(SHD) = 1 - 0.05 * SHD (0<SHD≤15), or 0.25 (SHD>15)
        SHD = SSH - SH
        f(SM) = 1 - exp(8.1 * (SM - 0.833))
        SH = 0.622 * e / (P - 0.378 * e)
        e = 6.112 * exp(17.67 * T / (T + 243.5))
    
        SHD: specific humidity deficit (g kg−1);
        SSH: saturated specific humidity;
        SH: specific humidity (g kg-1);
        SM: Soil moisture (unitless, 0-1);
        P:surface pressure (millibar, 1 millibar = 100 Pa);
        e: vapor pressure (millibar);
        T: Dew point temperature (°C), when calculate SSH, use air temperature (°C);
        */

        var SHD = Pmst.expression('0.622 * es / (P - 0.378 * es) - 0.622 * e / (P - 0.378 * e)',
            {
                'e': Tdew.expression('6.112 * exp(17.67 * T / (T + 243.5))', { 'T': Tdew }),
                'es': T.expression('6.112 * exp(17.67 * T / (T + 243.5))', { 'T': T }),
                'P': Pmst
            });
        var fw = SHD.expression('fSHD * fsm',
            {
                'fSHD': SHD.expression('1 - 0.05 * SHD', { 'SHD': SHD }).where(SHD.gt(15), 0.25),
                'fsm': sm.expression('1 - exp(8.1 * (SM - 0.833))', { 'SM': sm })
            });
    } else if (fw_param == '3-PG') {
        /*
        Ref: A generalised model of forest productivity using simplified concepts of radiation-use efficiency, carbon balance and partitioning.
    
        f(W) = min(f(SM), f(VPD))
        f(SM) = 1 / (1 + ((1 - (mr)/c) ** h)
        mr = SWC / SWCmax
        f(VPD) = exp(-k * VPD)
    
        mr: moisture ratio;
        SWC: soil water content (mm);
        SWCmax: maximum soil water content (mm);
        c, h: default to 0.55, 6;
        k: empirical coefficient, describing the relationship between stomatal and canopy conductance and D, default to 2.5;
        VPD: vapor pressure deficit (kPa);
        */

        var mr = sm.divide(sm_max);
        var fw = mr.expression('1 / (1 + ((1 - mr / 0.55) ** 6))', { 'mr': mr });
    } else if (fw_param == 'VPM') {
        // f(W) = (1 + LSWI) / (1 + LSWImax)
        if (interval <= 8 || end.difference(start, 'day').lte(8)) {
            var end8 = start.advance(8, 'day');
        }
        var nir_swir = ee.ImageCollection("MODIS/061/MOD09A1")
            .filter(ee.Filter.date(start, end8))
            .select(['sur_refl_b02', 'sur_refl_b06'])
            .map(allign_image)
            .mean();
        var lswi = nir_swir.expression('(nir - swir) / (nir + swir)', { 'nir': nir_swir.select('sur_refl_b02'), 'swir': nir_swir.select('sur_refl_b06') });
        var fw = lswi.expression('(1 + LSWI) / (1 + LSWImax)', { 'LSWI': lswi, 'LSWImax': lswi_max });
    } else if (fw_param == 'TEC') {
        // f(W) = E / Ep, E: actual evapotranspiration, Ep: potential evapotranspiration
        var fw = actual_eva.divide(potential_eva).clamp(0, 1);
    } else if (fw_param == 'CASA | MuSyQ') {
        // f(W) = 0.5 + 0.5 * E / Ep, E: actual evapotranspiration, Ep: potential evapotranspiration
        var fw = actual_eva.divide(potential_eva).multiply(0.5).add(0.5).clamp(0, 1);
    } else if (fw_param == 'TCF') {
        // f(W) = (SM - SMmin) / (SMmax - SMmin), range [0, 1]
        var fw = sm.expression('(SM - SMmin) / (SMmax - SMmin)', { 'SM': sm, 'SMmin': sm_min, 'SMmax': sm_max }).clamp(0, 1);
    } else if (fw_param == 'Wang-model') {
        // f(W) = (SWC - SWCmin) / (SWCmax - SWCmin), SWC: soil water content (mm)
        var fw = sm.expression('(SWC - SWCmin) / (SWCmax - SWCmin)', { 'SWC': sm, 'SWCmin': sm_min, 'SWCmax': sm_max }).clamp(0, 1);
    } else if (fw_param == 'CI-EF') {
        /*
        Ref: Improvements of the MODIS Gross Primary Productivity model based on a comprehensive uncertainty assessment over the Brazilian Amazonia.
    
        f(W) = (EF - EFmin) / (EFmax - EFmin)
        EF = LE / (H + LE)
    
        EF: plant evapotranspiration fraction. EFmin = 0.1, EFmax = 0.75
        LE: latent heat flux (W m-2)
        H: sensible heat flux (W m-2)
        
        */

        var le = ee.ImageCollection("NASA/GSFC/MERRA/lnd/2")
            .filter(ee.Filter.date(start, end))
            .select(["LHLAND"])
            .map(allign_image)
            .mean();
        var h = ee.ImageCollection("NASA/GSFC/MERRA/lnd/2")
            .filter(ee.Filter.date(start, end))
            .select(["SHLAND"])
            .map(allign_image)
            .mean();
        var ef = le.expression('LE / (H + LE)', { 'LE': le, 'H': h });
        var fw = ef.expression('(EF - EFmin) / (EFmax - EFmin)', { 'EF': ef, 'EFmin': 0.1, 'EFmax': 0.75 }).clamp(0, 1);
    } else if (fw_param == 'None') {
        var fw = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
    };
    return fw;
}

function calculate_fvpd_co2(fvpd_co2_param, T, Tk, vpd, vpd_min, vpd_max, Pmst, co2, resolution, lulc) {
    if (fvpd_co2_param == 'P-model') {
        /*
        f(VPD_CO2) = f(CO2) = mj * sqrt(1 - (0.41 / mj).power(2/3)) // P-model only got f(CO2), ignoring f(VPD)
        mj = (ci - gs) / (ci + 2 * gs)
        ci = chi * co2
    
        mj: Factor in the light-limited assimilation rate function (unitless);
        ci: Leaf-internal CO2 partial pressure (Pa);
        gs: photorespiratory CO2 compensation point in absence of dark respiration
        chi: pronounce of χ, Optimal ratio of leaf internal to ambient CO2 (unitless);
        co2: Ambient CO2 partial pressure (Pa);
        */

        var viscosity_of_water = calculateViscosityOfWater(T);
        var chi_gs = calculate_chi('pmodel', T, Tk, vpd, viscosity_of_water, Pmst, co2);
        var ci = chi_gs[0].multiply(co2);
        var mj = ci.expression('(ci - gs) / (ci + 2 * gs)', { 'ci': ci, 'gs': chi_gs[1] });
        var fco2 = mj.expression('mj * sqrt(1 - (0.41 / mj) ** (2/3))', { 'mj': mj });
        var fvpd_co2 = fco2;
    } else if (fvpd_co2_param == '3-PG') {
        // VPD's unit in '3-PG' is kPa (Ref. A generalised model of forest productivity using simplified concepts of radiation-use efficiency, carbon balance and partitioning. section 2.1. We used a value of kg = 2.5 (with vpd in KPa)).
        var fvpd_co2 = vpd.expression('exp(-0.25 * VPD)', { 'VPD': vpd });  // Change kg from original 2.5 to current 0.25, otherwise the value of fvpd_co2 will be too small (< 0.01).
    } else if (fvpd_co2_param == 'CFLUX | MOD17 | TL-LUE') {
        // f(VPD) = (VPDmax - VPD) / (VPDmax - VPDmin)
        var fvpd_co2 = vpd.expression('(VPDmax - VPD) / (VPDmax - VPDmin)', { 'VPD': vpd, 'VPDmax': vpd_max, 'VPDmin': vpd_min }).clamp(0, 1);
    } else if (fvpd_co2_param == 'C-Fix') {
        /*
        Ref: Estimation of carbon mass fluxes over Europe using the C-Fix model and Euroflux data, section 2.1.3 and A.2
    
        CO2fert = F(CO2) / F(CO2ref)
        F(CO2) = Vmax * ([CO2] - [O2] / (2 * τ)) / (Km * (1 + [O2] / K0) + [CO2])
        F(CO2ref) = Vmax * ([CO2ref] - [O2] / (2 * τ)) / (Km * (1 + [O2] / K0) + [CO2ref])
        τ = Aτ * exp(-Eaτ / Rg / T), Aτ = 7.87*1e-5, Eaτ = -42869.9
        Km = A * exp(-Ea / Rg / T)
        if T > 15 degC A = 2.419 * 1e13, Ea = 59.4 KJ mol-1
        else A = 1.976 * 1e22, Ea = 109.6 KJ mol-1
        K0 = A0 * exp(-Ea0 / Rg / T), A0 = 8240, Ea0 = 13913.5
    
        CO2fert: Normalised CO2 fertilisation factor, formalised according to Veroustraete (1994); no fertilisation means value equal to 1; fertilisation means values larger than 1, unitless;
        F(CO2): CO2 assimilation rate;
        [CO2]: CO2 concentration in the mesophyll tissue of leaves;
        [O2]: constant (assumed, not mentioned is paper), equal to 20.9;
        F(CO2ref): CO2 assimilation rate at reference conditions;
        [CO2ref]: CO2-concentration occurring in the reference year 1833, equal to 281 ppmv;
        τ: CO2/O2 specificity ratio, unitless;
        Km: Affinity constant for CO2 of Rubisco, unit is [% CO2];
        K0: Inhibition constant for O2, unit is [% O2];
        Vmax: Maximal photosynthetic rate, unit is g c m-2 s-1;
        Rg: Gas constant, equal to 8.31 J mol-1 K-1;
        T: Temperature, unit is K;
        */

        var k0 = Tk.expression('A0 * exp(-Ea0 / Rg / T)', { 'A0': 8240, 'Ea0': 13913.5, 'Rg': 8.31, 'T': Tk });
        var A = ee.Image(1.976 * 1e22).where(Tk.gt(288.15), 2.419 * 1e13);
        var Ea = ee.Image(109.6 * 1e3).where(Tk.gt(288.15), 59.4 * 1e3);
        var km = Tk.expression('A * exp(-Ea / Rg / T)', { 'A': A, 'Ea': Ea, 'Rg': 8.31, 'T': Tk });
        var tau = Tk.expression('A * exp(-Ea / Rg / T)', { 'A': 7.87 * 1e-5, 'Ea': -42869.9, 'Rg': 8.31, 'T': Tk });
        var Vmax = 1; // Vmax will be ignored by dividing F(CO2) and F(CO2ref)
        var Fco2 = tau.expression('Vmax * (CO2 - O2 / (2 * tau)) / (km * (1 + O2 / k0) + CO2)', { 'Vmax': Vmax, 'CO2': co2, 'O2': 20.9, 'tau': tau, 'km': km, 'k0': k0 });
        var Fco2ref = tau.expression('Vmax * (CO2ref - O2 / (2 * tau)) / (km * (1 + O2 / k0) + CO2ref)', { 'Vmax': Vmax, 'CO2ref': 281, 'O2': 20.9, 'tau': tau, 'km': km, 'k0': k0 });
        var fvpd_co2 = Fco2.divide(Fco2ref);
    } else if (fvpd_co2_param == 'EC-LUE-V3') {
        /*
        Ref: https://www.science.org/doi/full/10.1126/sciadv.aax1396
        In original EC-LUE-V3 model, f(CO2), f(T), f(VPD) are independent.
        Here, we combine f(CO2) and f(VPD) together as final f(VPD).
    
        f(VPD_CO2) = f(CO2) * min(f(T), f(VPD))
        f(CO2) = (Ci - theta) / (Ci + 2 * theta)
        f(VPD) = VPD0 / (VPD + VPD0)
        Ci = chi * co2
    
        Ci: Leaf-internal CO2 partial pressure (Pa);
        theta: pronuncing θ, the CO2 compensation point in the absence of dark respiration (ppm).
        chi: pronounce of χ, Optimal ratio of leaf internal to ambient CO2 (unitless);
        co2: atmospheric CO2 concentration (ppm);
        */

        var viscosity_of_water = calculateViscosityOfWater(T);
        var chi = calculate_chi('ec-lue', T, Tk, vpd, viscosity_of_water, Pmst, co2);
        var Ci = chi.multiply(co2);
        var theta_vpd0 = calculate_emax_theta_vpd0(lulc);
        var fco2 = Ci.expression('(Ci - theta) / (Ci + 2 * theta)', { 'Ci': Ci, 'theta': theta_vpd0[1] });
        var fvpd = vpd.expression('VPD0 / (VPD + VPD0)', { 'VPD0': theta_vpd0[2], 'VPD': vpd });
        // var fvpd_co2 = fco2.multiply(fvpd.min(ft));
        var fvpd_co2 = fco2.multiply(fvpd);
    } else if (fvpd_co2_param == 'CCW') {
        // f(VPD) = exp(-K2 * (VPD - VPDmin)), unit of VPD is hPa
        var k2 = calculate_bplut_ccw(lulc)[3];
        var fvpd_co2 = vpd.expression('exp(-K2 * (VPD - VPDmin))', { 'K2': k2, 'VPD': vpd.multiply(10), 'VPDmin': vpd_min.multiply(10) });
    } else if (fvpd_co2_param == 'None') {
        var fvpd_co2 = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
    };
    return fvpd_co2;
}

// calculate VPD
function calculate_vpd(Pmst, T, Tdew) {
    /*
    Formula from https://www.science.org/doi/full/10.1126/sciadv.aax1396.
    SVP and AVP are saturated vapor pressure and actual vapor pressure (hPa), respectively. Tdew is the dew point temperature (°C). RH is the land relative humidity (%).
    T is the land air temperature (°C). Z is the altitude (m). Pmst is the air pressure (hPa), and Pmsl is the air pressure at mean sea level (1013.25 hPa).
    Attention: The unit here for Pmst, svp and vpd are hPa, We need to convert it to kPa eventually.
    */

    var func_w = Pmst.expression('1 + 7 * 1e-4 + 3.46 * 1e-6 * Pmst', { 'Pmst': Pmst }); // Accoording to this formula, the Pmst (range from 800-1200 hPa) has almost no effect on the func_w and svp???
    var svp = func_w.expression('6.112 * func_w * exp(17.76 * T / (T + 243.5))', { 'func_w': func_w, 'T': T });
    var rh = T.expression('exp(Tdew * 17.629 / (Tdew + 273.3) - (T * 17.629 / (T + 237.3)))', { 'Tdew': Tdew, 'T': T });
    var vpd = svp.expression('svp * (1 - rh)', { 'svp': svp, 'rh': rh });
    return vpd.divide(10); // hPa to kPa
}

function calculate_chi(model, T, Tk, vpd, viscosity_of_water, Pmst, co2) {
    /*
    # Description: calculate χ (chi), the optimal ratio of leaf internal to ambient CO2.
    # Reference: rpmodel package manual (https://cran.r-project.org/web/packages/rpmodel/rpmodel.pdf). Collected from function ftemp_arrh, gammastar, kmm, rpmodel
    Formular:
    if model == "ec-lue":
        χ = Xi / (Xi + sqrt(vpd))
        Xi = np.sqrt(356.51 * K / 1.6 / η∗)
    elif model == "pmodel":
        χ = gs / ca + (1 - gs / ca) *  Xi / (Xi + sqrt(vpd))
        gs = gs0 * f(T, ΔHa) * pz / 101325 = 4.332 * exp(37830 / 8.314 * (1 / Tk - 1 / 298.15))
        Xi = sqrt(146 * (K + gs) / 1.6 / η∗)
    η∗ = Eta(Tk) / Eta(25degC)
    K = Kc * (1 + pO2 / Ko)
    pO2 = 0.209476 * Pmst
    Kc = Kc25 * f(T, ΔHkc) = Kc25 * exp(∆Hkc / R * (1 / Tk - 1 / Tkref)) = 39.97 * exp(79430 / 8.314 * (1 / Tk - 1 / 298.15) = 39.97 * exp(9553.76 / Tk - 32.04)
    Ko = K025 * f(T, ΔHko) = 27480 * exp(4375.75 / Tk - 14.68)
    f = exp(ΔH/R * (1/Tk - 1/Tkref))
  
    χ: Optimal ratio of leaf internal to ambient CO2 (unitless);
    gs: photorespiratory CO2 compensation point in absence of dark respiration;
    gs0: gs value at standard temperature (T0 = 25deg C) and atmospheric pressure (p0 = 101325Pa), constant;
    Kc: the Michaelis-Menten constant for CO2 (Pa);
    Ko: the Michaelis-Menten constant for O2 (Pa);
    f: the temperature scaling function;
    R: the universal gas constant (J mol-1 K-1);
    Tk: the air temperature (K);
    Tkref: the reference temperature (K);
    pO2: the partial pressure of oxygen (Pa), calculated as 0.209476 * Atmospheric pressure (Pa);
    Pmst: the atmospheric pressure (Pa);
    K: Michaelis Menten coefficient of Rubisco-limited assimilation as a function of tem-perature and atmospheric pressure;
    η∗: Change in the viscosity of water, relative to its value at 25 deg C (unitless);
    vpd: the vapour pressure deficit (Pa);
  
    Constants: gs0(4.332 Pa), ΔHa(37830 J mol-1), ∆Hkc(79430 J mol-1), ∆Hko(36380 J mol-1), Kc25(39.97 Pa), and Ko25(27480 Pa), R(8.314 J mol-1 K-1), Tkref(298.15 K), p0(101325 Pa)
    */

    var vpd = vpd.multiply(1000); // kPa to Pa
    var Kc = Tk.expression('39.97 * exp(9553.76 / Tk - 32.04)', { 'Tk': Tk });
    var Ko = Tk.expression('27480 * exp(4375.75 / Tk - 14.68)', { 'Tk': Tk });
    var pO2 = Pmst.multiply(0.209476);
    var K = Kc.expression('Kc * (1 + pO2 / Ko)', { 'Kc': Kc, 'pO2': pO2, 'Ko': Ko });
    var viscosity_of_water = calculateViscosityOfWater(T);
    if (model == "ec-lue") {
        var Xi = K.expression('sqrt(356.51 * K / 1.6 / Eta)', { 'K': K, 'Eta': viscosity_of_water });
        var χ = Xi.expression('Xi / (Xi + sqrt(vpd))', { 'Xi': Xi, 'vpd': vpd });
        return χ;
    } else if (model == "pmodel") {
        var gs = Tk.expression('4.332 * exp(37830 / 8.314 * (1 / Tk - 1 / 298.15))', { 'Tk': Tk });
        var Xi = K.expression('sqrt(91.25 * (K + gs) / Eta)', { 'K': K, 'gs': gs, 'Eta': viscosity_of_water });
        // let a = gs / ca, b = Xi / (Xi + sqrt(vpd)), χ = a + b - a * b
        var a = gs.divide(co2);
        var b = Xi.expression('Xi / (Xi + sqrt(vpd))', { 'Xi': Xi, 'vpd': vpd });
        var χ = a.expression('a + b - a * b', { 'a': a, 'b': b });
        return [χ, gs];
    }
}

function calculateViscosityOfWater(T) {
    // https://byjus.com/chemistry/viscosity-of-water/
    var temperature = [50, 45, 40, 35, 30, 25, 20, 15, 10, 5];
    var viscosity = [0.5465, 0.5958, 0.6527, 0.7191, 0.7972, 0.89, 1.0016, 1.1375, 1.3059, 1.5182]; // unit: mPa.s
    var viscosity_of_water = ee.Image(0);
    for (var i = 0; i < temperature.length; i++) {
        viscosity_of_water = viscosity_of_water.where(T.lt(temperature[i]), ee.Image.constant(viscosity[i]));
    }
    return viscosity_of_water.divide(0.89); // return change in the viscosity of water, relative to its value at 25°C (unitless)
}

function calculate_emax_theta_vpd0(lulc) {
    var emax = ee.Image(0);
    var theta = ee.Image(0);
    var vpd0 = ee.Image(0);
    var value_dict = [{ "class_value": 1, "Vegetation_Type": "ENF", "emax": 3.16, "theta": 20.69, "vpd0": 0.69 },
    { "class_value": 2, "Vegetation_Type": "EBF", "emax": 3.92, "theta": 30.0, "vpd0": 0.29 },
    { "class_value": 3, "Vegetation_Type": "DBF", "emax": 2.02, "theta": 60.0, "vpd0": 0.54 },
    { "class_value": 4, "Vegetation_Type": "DBF", "emax": 2.02, "theta": 60.0, "vpd0": 0.54 },
    { "class_value": 5, "Vegetation_Type": "MF", "emax": 3.46, "theta": 52.0, "vpd0": 0.41 },
    { "class_value": 6, "Vegetation_Type": "SHR", "emax": 0.75, "theta": 37.0, "vpd0": 1.13 },
    { "class_value": 7, "Vegetation_Type": "SHR", "emax": 0.75, "theta": 37.0, "vpd0": 1.13 },
    { "class_value": 8, "Vegetation_Type": "SAV", "emax": 1.54, "theta": 41.0, "vpd0": 1.04 },
    { "class_value": 9, "Vegetation_Type": "SAV", "emax": 1.54, "theta": 41.0, "vpd0": 1.04 },
    { "class_value": 10, "Vegetation_Type": "GRA", "emax": 3.32, "theta": 75.0, "vpd0": 1.21 },
    { "class_value": 11, "Vegetation_Type": "WET", "emax": 2.65, "theta": 68.0, "vpd0": 1.16 },
    { "class_value": 12, "Vegetation_Type": "CRO", "emax": 2.85, "theta": 64, "vpd0": 0.89 },
    { "class_value": 13, "Vegetation_Type": "Urban", "emax": 0, "theta": 0, "vpd0": 0 },
    { "class_value": 14, "Vegetation_Type": "CRO", "emax": 2.85, "theta": 64, "vpd0": 0.89 },
    { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "emax": 0, "theta": 0, "vpd0": 0 },
    { "class_value": 16, "Vegetation_Type": "Barren", "emax": 0, "theta": 0, "vpd0": 0 },
    { "class_value": 17, "Vegetation_Type": "Water Bodies", "emax": 0, "theta": 0, "vpd0": 0 }];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var emax = emax.where(lulc.eq(value["class_value"]), value["emax"]);
        var theta = theta.where(lulc.eq(value["class_value"]), value["theta"]);
        var vpd0 = vpd0.where(lulc.eq(value["class_value"]), value["vpd0"]);
    }
    return [emax, theta, vpd0];
}

function calculate_bplut_modis(lulc) {
    var Tmin_min = ee.Image(0);
    var Tmin_max = ee.Image(0);
    var lue_max = ee.Image(0);
    var value_dict = [
        { "class_value": 1, "Vegetation_Type": "ENF", "Tmin_min": -8, "Tmin_max": 8.31, "lue_max": 0.962 },
        { "class_value": 2, "Vegetation_Type": "EBF", "Tmin_min": -8, "Tmin_max": 9.09, "lue_max": 1.268 },
        { "class_value": 3, "Vegetation_Type": "DNF", "Tmin_min": -8, "Tmin_max": 9.09, "lue_max": 1.086 },
        { "class_value": 4, "Vegetation_Type": "DBF", "Tmin_min": -6, "Tmin_max": 9.94, "lue_max": 1.165 },
        { "class_value": 5, "Vegetation_Type": "MF", "Tmin_min": -8, "Tmin_max": 9.50, "lue_max": 1.051 },
        { "class_value": 6, "Vegetation_Type": "CSH", "Tmin_min": -8, "Tmin_max": 8.61, "lue_max": 1.281 },
        { "class_value": 7, "Vegetation_Type": "OSH", "Tmin_min": -8, "Tmin_max": 8.80, "lue_max": 0.841 },
        { "class_value": 8, "Vegetation_Type": "SAV", "Tmin_min": -8, "Tmin_max": 11.39, "lue_max": 1.239 },
        { "class_value": 9, "Vegetation_Type": "SAV", "Tmin_min": -8, "Tmin_max": 11.39, "lue_max": 1.206 },
        { "class_value": 10, "Vegetation_Type": "GRA", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 0.860 },
        { "class_value": 11, "Vegetation_Type": "WET", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 0.860 },
        { "class_value": 12, "Vegetation_Type": "CRO", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 1.044 },
        { "class_value": 13, "Vegetation_Type": "Urban", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0 },
        { "class_value": 14, "Vegetation_Type": "CRO", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 1.044 },
        { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0 },
        { "class_value": 16, "Vegetation_Type": "Barren", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0 },
        { "class_value": 17, "Vegetation_Type": "Water Bodies", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0 }
    ];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var Tmin_min = Tmin_min.where(lulc.eq(value["class_value"]), value["Tmin_min"]);
        var Tmin_max = Tmin_max.where(lulc.eq(value["class_value"]), value["Tmin_max"]);
        var lue_max = lue_max.where(lulc.eq(value["class_value"]), value["lue_max"]);
    };
    return [Tmin_min, Tmin_max, lue_max];
}

function calculate_bplut_CI_EF(lulc) {
    var Tmin_min = ee.Image(0);
    var Tmin_max = ee.Image(0);
    var lue_max = ee.Image(0);
    var lue_cs = ee.Image(0);
    var value_dict = [
        { "class_value": 1, "Vegetation_Type": "ENF", "Tmin_min": -8, "Tmin_max": 8.31, "lue_max": 0.96, "lue_cs": 0.78 },
        { "class_value": 2, "Vegetation_Type": "EBF", "Tmin_min": -8, "Tmin_max": 9.09, "lue_max": 2.23, "lue_cs": 0.78 },
        { "class_value": 3, "Vegetation_Type": "DNF", "Tmin_min": -8, "Tmin_max": 9.09, "lue_max": 2.23, "lue_cs": 0.78 },
        { "class_value": 4, "Vegetation_Type": "DBF", "Tmin_min": -6, "Tmin_max": 9.94, "lue_max": 2.23, "lue_cs": 0.78 },
        { "class_value": 5, "Vegetation_Type": "MF", "Tmin_min": -8, "Tmin_max": 9.50, "lue_max": 1.40, "lue_cs": 0.65 },
        { "class_value": 6, "Vegetation_Type": "CSH", "Tmin_min": -8, "Tmin_max": 8.61, "lue_max": 1.28, "lue_cs": 1.28 }, // The papaer did not mention lue_cs of Shurbland, so set it equal to lue_max
        { "class_value": 7, "Vegetation_Type": "OSH", "Tmin_min": -8, "Tmin_max": 8.80, "lue_max": 0.84, "lue_cs": 0.84 },
        { "class_value": 8, "Vegetation_Type": "SAV", "Tmin_min": -8, "Tmin_max": 11.39, "lue_max": 0.80, "lue_cs": 0.48 },
        { "class_value": 9, "Vegetation_Type": "SAV", "Tmin_min": -8, "Tmin_max": 11.39, "lue_max": 0.80, "lue_cs": 0.48 },
        { "class_value": 10, "Vegetation_Type": "GRA", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 1.04, "lue_cs": 0.55 },
        { "class_value": 11, "Vegetation_Type": "WET", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 1.04, "lue_cs": 0.55 },
        { "class_value": 12, "Vegetation_Type": "CRO", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 0.73, "lue_cs": 0.37 },
        { "class_value": 13, "Vegetation_Type": "Urban", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0, "lue_cs": 0 },
        { "class_value": 14, "Vegetation_Type": "CRO", "Tmin_min": -8, "Tmin_max": 12.02, "lue_max": 0.73, "lue_cs": 0.37 },
        { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0, "lue_cs": 0 },
        { "class_value": 16, "Vegetation_Type": "Barren", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0, "lue_cs": 0 },
        { "class_value": 17, "Vegetation_Type": "Water Bodies", "Tmin_min": 0, "Tmin_max": 0, "lue_max": 0, "lue_cs": 0 }
    ];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var Tmin_min = Tmin_min.where(lulc.eq(value["class_value"]), value["Tmin_min"]);
        var Tmin_max = Tmin_max.where(lulc.eq(value["class_value"]), value["Tmin_max"]);
        var lue_max = lue_max.where(lulc.eq(value["class_value"]), value["lue_max"]);
        var lue_cs = lue_cs.where(lulc.eq(value["class_value"]), value["lue_cs"]);
    };
    return [Tmin_min, Tmin_max, lue_max, lue_cs];
}

function calculate_bplut_MuSyQ(lulc) {
    var lue_su = ee.Image(0);
    var lue_sh = ee.Image(0);
    var Topt = ee.Image(0);
    var value_dict = [
        { "class_value": 1, "Vegetation_Type": "ENF", "lue_su": 0.678, "lue_sh": 3.020, "Topt": 15 },
        { "class_value": 2, "Vegetation_Type": "EBF", "lue_su": 0.706, "lue_sh": 3.079, "Topt": 25 },
        { "class_value": 3, "Vegetation_Type": "DNF", "lue_su": 0.403, "lue_sh": 1.700, "Topt": 15 },
        { "class_value": 4, "Vegetation_Type": "DBF", "lue_su": 0.680, "lue_sh": 3.030, "Topt": 20 },
        { "class_value": 5, "Vegetation_Type": "MF", "lue_su": 0.795, "lue_sh": 2.917, "Topt": 17 },
        { "class_value": 6, "Vegetation_Type": "CSH", "lue_su": 0.552, "lue_sh": 2.446, "Topt": 20 },
        { "class_value": 7, "Vegetation_Type": "OSH", "lue_su": 0.366, "lue_sh": 1.807, "Topt": 16 },
        { "class_value": 8, "Vegetation_Type": "SAV", "lue_su": 0.615, "lue_sh": 2.914, "Topt": 20 },
        { "class_value": 9, "Vegetation_Type": "SAV", "lue_su": 0.615, "lue_sh": 2.914, "Topt": 20 },
        { "class_value": 10, "Vegetation_Type": "GRA", "lue_su": 0.603, "lue_sh": 2.833, "Topt": 18 },
        { "class_value": 11, "Vegetation_Type": "WET", "lue_su": 0.603, "lue_sh": 2.914, "Topt": 18 },
        { "class_value": 12, "Vegetation_Type": "CRO", "lue_su": 1.114, "lue_sh": 2.913, "Topt": 26 },
        { "class_value": 13, "Vegetation_Type": "Urban", "lue_su": 0, "lue_sh": 0, "Topt": 0 },
        { "class_value": 14, "Vegetation_Type": "CRO", "lue_su": 1.114, "lue_sh": 2.913, "Topt": 26 },
        { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "lue_su": 0, "lue_sh": 0, "Topt": 0 },
        { "class_value": 16, "Vegetation_Type": "Barren", "lue_su": 0, "lue_sh": 0, "Topt": 0 },
        { "class_value": 17, "Vegetation_Type": "Water Bodies", "lue_su": 0, "lue_sh": 0, "Topt": 0 }
    ];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var lue_su = lue_su.where(lulc.eq(value["class_value"]), value["lue_su"]);
        var lue_sh = lue_sh.where(lulc.eq(value["class_value"]), value["lue_sh"]);
        var Topt = Topt.where(lulc.eq(value["class_value"]), value["Topt"]);
    };
    return [lue_su, lue_sh, Topt];
}

function calculate_bplut_ccw(lulc) {
    var lue_opt = ee.Image(0);
    var Topt = ee.Image(0);
    var K1 = ee.Image(0);
    var K2 = ee.Image(0);
    var value_dict = [
        { "class_value": 1, "Vegetation_Type": "ENF", "lue_opt": 2.516, "Topt": 27.8, "K1": 0.575, "K2": 0.091 },
        { "class_value": 2, "Vegetation_Type": "EBF", "lue_opt": 2.892, "Topt": 28.0, "K1": 0.745, "K2": 0.05 },
        { "class_value": 3, "Vegetation_Type": "DNF", "lue_opt": 2.233, "Topt": 28.0, "K1": 0.550, "K2": 0.032 },
        { "class_value": 4, "Vegetation_Type": "DBF", "lue_opt": 2.233, "Topt": 28.0, "K1": 0.550, "K2": 0.032 },
        { "class_value": 5, "Vegetation_Type": "MF", "lue_opt": 2.758, "Topt": 27.7, "K1": 0.755, "K2": 0.053 },
        { "class_value": 6, "Vegetation_Type": "CSH", "lue_opt": 2.209, "Topt": 26.7, "K1": 0.643, "K2": 0.095 },
        { "class_value": 7, "Vegetation_Type": "OSH", "lue_opt": 2.144, "Topt": 27.9, "K1": 0.789, "K2": 0.057 },
        { "class_value": 8, "Vegetation_Type": "SAV", "lue_opt": 2.516, "Topt": 27.8, "K1": 0.575, "K2": 0.057 },
        { "class_value": 9, "Vegetation_Type": "SAV", "lue_opt": 2.190, "Topt": 26.7, "K1": 0.643, "K2": 0.032 },
        { "class_value": 10, "Vegetation_Type": "GRA", "lue_opt": 2.397, "Topt": 16.9, "K1": 0.665, "K2": 0.088 },
        { "class_value": 11, "Vegetation_Type": "WET", "lue_opt": 0, "Topt": 0, "K1": 0, "K2": 0 },
        { "class_value": 12, "Vegetation_Type": "CRO", "lue_opt": 2.319, "Topt": 28.0, "K1": 0.254, "K2": 0.096 },
        { "class_value": 13, "Vegetation_Type": "Urban", "lue_opt": 0, "Topt": 0, "K1": 0, "K2": 0 },
        { "class_value": 14, "Vegetation_Type": "CRO", "lue_opt": 2.319, "Topt": 28.0, "K1": 0.254, "K2": 0.096 },
        { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "lue_opt": 0, "Topt": 0, "K1": 0, "K2": 0 },
        { "class_value": 16, "Vegetation_Type": "Barren", "lue_opt": 0, "Topt": 0, "K1": 0, "K2": 0 },
        { "class_value": 17, "Vegetation_Type": "Water Bodies", "lue_opt": 0, "Topt": 0, "K1": 0, "K2": 0 }
    ];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var lue_opt = lue_opt.where(lulc.eq(value["class_value"]), value["lue_opt"]);
        var Topt = Topt.where(lulc.eq(value["class_value"]), value["Topt"]);
        var K1 = K1.where(lulc.eq(value["class_value"]), value["K1"]);
        var K2 = K2.where(lulc.eq(value["class_value"]), value["K2"]);
    };
    return [lue_opt, Topt, K1, K2];
}

function generate_c3c4_map(lulc) {
    /*
    The vast majority of forests and shrubs can be considered as C3 plants [1]. For grasslands, we can determine whether they are C3 or C4 plants based on annual accumulated temperature and annual precipitation [2]. For crops, due to the lack of accurate and continuously updated global crop-distribution datasets [3], they are temporarily classified as a mixture of 50% C3 and 50% C4.
  
    The global distribution of C3 and C4 grasslands is influenced by various factors such as climate, temperature, precipitation, and soil conditions. Here are some general patterns of the global distribution of C3 and C4 grasslands:
    C3 Grasslands:
      C3 grasslands are typically found in cooler climates, including temperate regions and high latitudes.
      They are common in regions with higher annual precipitation and more favorable moisture conditions.
      Examples of C3 grasslands include the prairies of North America, the steppes of Eurasia, and the Pampas of South America.
    C4 Grasslands:
      C4 grasslands are predominantly found in warmer climates, including tropical and subtropical regions.
      They are often associated with drier environments, lower annual precipitation, and higher temperatures.
      Examples of C4 grasslands include the savannas of Africa, the cerrados of South America, and the grasslands of Australia.
    It's important to note that the distribution of C3 and C4 grasslands is not strictly defined by latitude or region. There can be transitional areas where both types of grasslands coexist or where the dominance of one type shifts gradually into the other. Additionally, factors like fire regimes, grazing pressure, and land-use practices can also influence the composition and distribution of C3 and C4 grasslands within a specific region.
  
    Look-up table for grassland [2]:
      Preliminary grassland-type PFT | Subdivided grassland-type PFT | Climate rules
      Grassland | C3 grass, arctic | GDD<400
      Grassland | C3 grass | GDD ≥ 400 and (Tw ≤ 22°C or six months Pmon ≤ 25 mm and for month Tmon > 22°C) (PS: Does 'six months Pmon ≤ 25 mm and for month Tmon > 22°C' means 'for Tmon > 22°C, the Pmon of the month is less than 25 mm, and this kind of month should be more than six in a year'?)
      Grassland | C4 grass | GDD ≥ 400 and Tc ≥ 22°C and driest month Pmon > 25 mm
      Grassland | Mixed C3/C4 grass | Other grasslands that do not meet the above rules
    
    GDD: growing-degree day. The annual accumulated temperature of days over 5°C. For example, if the temperature of a day is 10°C, then the GDD of this day is 5°C.
    Tw: The average temperature in the warmest month.
    Tc: The average temperature in the coldest month.
    Tmon: The monthly temperature.
    Pmon: The monthly precipitation.
  
    Ref:
      1. Developing a diagnostic model for estimating terrestrial vegetation gross primary productivity using the photosynthetic quantum yield and Earth Observation data.
      2. Global land projection based on plant functional types with a 1-km resolution under socio-climatic scenarios.
      3. Global distribution of C3 and C4 vegetation: Carbon cycle implications.
    */

    // 10: C3, 20: C4, 30: Mixed C3/C4, 40: Others
    var c3c4_map = ee.Image(40);
    var year = ee.Date(start_date_select.getValue()[0]).get('year');
    var gdd = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1))
        .select('temperature_2m')
        .reduce(ee.Reducer.sum())
        .subtract(278.15 * 365);  // Calculate GDD
    var dataset = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1));
    var Tmon = dataset.select('temperature_2m');
    var Pmon = dataset.select('total_precipitation_sum');
    var Tw = Tmon.max();
    var Tc = Tmon.min();
    var Pdry = Pmon.min();
    var Pmon_lte25_and_Tmon_gt22 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
        .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1))
        .map(function (image) {
            var Tmon = image.select('temperature_2m');
            var Pmon = image.select('total_precipitation_sum');
            var image = image.updateMask(Pmon.lte(0.025)).updateMask(Tmon.gt(273.15 + 22));
            return image.multiply(0).add(1);  // Eligible months got a value of 1.
        })
        .select(1)
        .sum();
    var condition0 = gdd.lt(400);
    var condition1 = gdd.gte(400).and(Tw.lte(273.15 + 22).or(Pmon_lte25_and_Tmon_gt22.gte(6)));
    var condition2 = gdd.gte(400).and(Tc.gte(273.15 + 22)).or(Pdry.gt(0.015));  // If we follow the definition in [2], there is almost no C4 grassland in the calculation result. Therefore, the condition has been modified here. "and(Pdry.gt(0.025))" has been changed to "or(Pdry.gt(0.015))".
    var c3c4_map = c3c4_map.where(condition0, 10);
    var c3c4_map = c3c4_map.where(condition1, 10);
    var c3c4_map = c3c4_map.where(condition2, 20);
    var c3c4_map = c3c4_map.where(c3c4_map.neq(10).or(c3c4_map.neq(20)), 30);

    var value_dict = [
        { "class_value": 1, "Vegetation_Type": "ENF", "c3c4": 10 },
        { "class_value": 2, "Vegetation_Type": "EBF", "c3c4": 10 },
        { "class_value": 3, "Vegetation_Type": "DNF", "c3c4": 10 },
        { "class_value": 4, "Vegetation_Type": "DBF", "c3c4": 10 },
        { "class_value": 5, "Vegetation_Type": "MF", "c3c4": 10 },
        { "class_value": 6, "Vegetation_Type": "CSH", "c3c4": 10 },
        { "class_value": 7, "Vegetation_Type": "OSH", "c3c4": 10 },
        { "class_value": 8, "Vegetation_Type": "SAV", "c3c4": 10 },
        { "class_value": 9, "Vegetation_Type": "SAV", "c3c4": 10 },
        // { "class_value": 10, "Vegetation_Type": "GRA", "c3c4": 20 },
        { "class_value": 11, "Vegetation_Type": "WET", "c3c4": 40 },
        { "class_value": 12, "Vegetation_Type": "CRO", "c3c4": 30 },
        { "class_value": 13, "Vegetation_Type": "Urban", "c3c4": 40 },
        { "class_value": 14, "Vegetation_Type": "CRO", "c3c4": 30 },
        { "class_value": 15, "Vegetation_Type": "Permanent Snow and Ice", "c3c4": 40 },
        { "class_value": 16, "Vegetation_Type": "Barren", "c3c4": 40 },
        { "class_value": 17, "Vegetation_Type": "Water Bodies", "c3c4": 40 }
    ];
    for (var i = 0; i < value_dict.length; i++) {
        var value = value_dict[i];
        var c3c4_map = c3c4_map.where(lulc.eq(value["class_value"]), value["c3c4"]);
    };
    return c3c4_map;
}

function allign_image(image) {
    if (image.projection().crs() != 'EPSG:4326') {
        image = image.reproject('EPSG:4326');
    }
    image = image.clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
    return image;
}

function calculate_ci(image) {
    var result = image.expression('1 - actural_radiance / potential_radiance', { 'actural_radiance': image.select('SWGNT'), 'potential_radiance': image.select('SWTNT') });
    return result;
};

function calculate_lswi(image) {
    var result = image.expression('(nir - swir) / (nir + swir)', { 'nir': image.select('sur_refl_b02'), 'swir': image.select('sur_refl_b06') });
    return result;
};

function calculate_invariant_variables(start_date, fpar_param, ft_param, fw_param, fvpd_co2_param, lue_ci_param) {
    // Variables that are not related to time
    var year = start_date.get('year');
    var long_period_start = ee.Date.fromYMD(year.add(-2), 1, 1);
    var long_period_end = ee.Date.fromYMD(year.add(2), 1, 1);
    var co2_mm_gl = ee.FeatureCollection('projects/ee-849240929chen/assets/co2_mm_gl');
    if (fpar_param == 'CASA') {
        // Calculate SRmax and SRmin for each pixel, once for all. SR = (1 + NDVI) / (1 - NDVI)
        var sr = ee.ImageCollection("MODIS/MOD09GA_006_NDVI")
            .filterDate(long_period_start, long_period_end)
            .select('NDVI')
            .map(allign_image)
            .map(function (image) {
                return image.expression('(1 + ndvi) / (1 - ndvi)', { 'ndvi': image });
            });
        var sr_max = sr.max().clamp(4, 7);
        var sr_min = sr.min().clamp(0.5, 1.5);
    }
    if (lue_ci_param == 'CFLUX | CI-EF' || lue_ci_param == 'Wang-model') {
        // calculate maximum and minimum CI for each pixel, once for all.
        var ci_long_period = ee.ImageCollection('NASA/GSFC/MERRA/rad/2')
            .filterDate(long_period_start, long_period_end)
            .map(allign_image)
            .map(calculate_ci);
        var ci_max = ci_long_period.max();
        var ci_min = ci_long_period.min();
    }
    if (ft_param == 'CASA') {
        // Calculate Topt for each pixel, once for all.
        var ndvi = ee.ImageCollection("MODIS/061/MOD13A1")
            .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1))
            .select('NDVI');
        var ndvi_max = ndvi.max();
        var ndvi_max_month = ndvi.map(function (image) {
            // The best way to create a zero-value copy of an image is to use image.subtract(image). Using ee.Image(0) will not copy metadata such as coordinate systems, shape, and data type from the original image.
            return image.subtract(image).where(image.eq(ndvi_max), image.date().get('month'));
        }).map(allign_image)
            .max();
        var Topt_CASA = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
            .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 1, 1))
            .select('temperature_2m')
            .map(allign_image)
            .map(function (image) {
                return image.subtract(image).where(ndvi_max_month.eq(image.date().get('month')), image)
            }).max()
            .subtract(273.15);
    }
    if (fw_param == 'VPM') {
        // lswi_max: the maximum LSWI during the growing season for each pixel
        var lswi_max = ee.ImageCollection("MODIS/061/MOD09A1")
            .filterDate(long_period_start, long_period_end)
            .map(allign_image)
            .map(calculate_lswi)
            .max();
    } else if (fw_param == 'TCF' || fw_param == '3-PG' || fw_param == 'Wang-model') {
        // calculate maximum and minimum soil moisture for each pixel, once for all.
        var average = function (image) {
            return image.expression('(SM10 + SM40) / 2', { 'SM10': image.select('SoilMoi00_10cm_tavg'), 'SM40': image.select('SoilMoi10_40cm_tavg') });
        }
        var sm_long_peroid = ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001")
            .filterDate(long_period_start, long_period_end)
            .map(allign_image)
            .map(average);
        var sm_max = sm_long_peroid.max();
        var sm_min = sm_long_peroid.min();
    }

    if (fvpd_co2_param == 'CFLUX | MOD17 | TL-LUE' | fvpd_co2_param == 'CCW') {
        // calculate maximum and minimum VPD for each pixel, once for all.
        /* By TerraClimate, the resulting vpd_max and vpd_min are not consistent with vpd calculated by function calculate_vpd.
        var vpd_long_period = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
          .filterDate(long_period_start, long_period_end)
          .select('vpd')
          .map(allign_image);
        var vpd_max = vpd_long_period.max().multiply(1e-2); // scale factor: 0.01, unit: kPa
        var vpd_min = vpd_long_period.min().multiply(1e-2);
        */

        // By ERA5-Land
        function calculate_vpd_long_peroid(image) {
            var T = image.select('temperature_2m').subtract(273.15);
            var Tdew = image.select('dewpoint_temperature_2m').subtract(273.15);
            var svp = T.expression('6.112 * exp(17.76 * T / (T + 243.5))', { 'T': T });
            var rh = T.expression('exp(Tdew * 17.629 / (Tdew + 273.3) - (T * 17.629 / (T + 237.3)))', { 'Tdew': Tdew, 'T': T });
            var vpd = svp.expression('svp * (1 - rh)', { 'svp': svp, 'rh': rh });
            return vpd.divide(10); // hPa to kPa
        }

        var vpd_long_peroid = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
            .filterDate(long_period_start, long_period_end)
            .map(allign_image)
            .map(calculate_vpd_long_peroid);
        var vpd_max = vpd_long_peroid.max();
        var vpd_min = vpd_long_peroid.min();
    }
    return [co2_mm_gl, lswi_max, sm_max, sm_min, vpd_max, vpd_min, ci_max, ci_min, Topt_CASA, sr_max, sr_min];
}

function calculate_shared_variables(start, end, co2_mm_gl) {
    // Varaibles that are reused in different functions
    var Tdew = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filter(ee.Filter.date(start, end))
        .select(['dewpoint_temperature_2m'])
        .map(allign_image)
        .mean()
        .subtract(273.15); // convert K to °C
    var Tk = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filter(ee.Filter.date(start, end))
        .select(['temperature_2m'])
        .map(allign_image)
        .mean();
    var T = Tk.subtract(273.15); // convert K to °C
    var Pmst = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filter(ee.Filter.date(start, end))
        .select(['surface_pressure'])
        .map(allign_image)
        .mean()
        .multiply(1e-2); // convert Pa to hPa
    var sm = ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001")
        .filter(ee.Filter.date(start.update(null, null, 1), start.update(null, null, 2)))
        .select(['SoilMoi00_10cm_tavg', 'SoilMoi10_40cm_tavg'])
        .map(allign_image)
        .mean()
        .reduce(ee.Reducer.mean());
    var potential_eva = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filter(ee.Filter.date(start, end))
        .select(['potential_evaporation_sum'])
        .map(allign_image)
        .mean()
        .multiply(0.05);  // too much compared to actual_eva, so multiply 0.05 to reduce the value
    var actual_eva = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
        .filter(ee.Filter.date(start, end))
        .select(['evaporation_from_vegetation_transpiration_sum'])
        .map(allign_image)
        .mean();
    var vpd = calculate_vpd(Pmst, T, Tdew); // unit: kPa
    var co2 = ee.Image.constant(co2_mm_gl.filterMetadata('year', 'equals', start.get('year')).filterMetadata('month', 'equals', start.get('month')).first().toDictionary().get('average'));
    var lulc = ee.ImageCollection('MODIS/061/MCD12Q1')
        .filterDate(start.advance(-1, 'year'), start)
        .select('LC_Type1')
        .map(allign_image)
        .first();
    var c3c4_map = generate_c3c4_map(lulc);
    var ci = ee.ImageCollection('NASA/GSFC/MERRA/rad/2')
        .filter(ee.Filter.date(start, end))
        .map(allign_image)
        .map(calculate_ci)
        .mean();
    return [T, Tk, Tdew, Pmst, sm, potential_eva, actual_eva, vpd, co2, lulc, c3c4_map, ci];
}


function main(start_date, end_date, resolution, interval, par_param, fpar_param, lue_ci_param, ft_param, fw_param, fvpd_co2_param) {
    var iv = calculate_invariant_variables(start_date, fpar_param, ft_param, fw_param, fvpd_co2_param, lue_ci_param);
    var co2_mm_gl = iv[0], lswi_max = iv[1], sm_max = iv[2], sm_min = iv[3], vpd_max = iv[4], vpd_min = iv[5], ci_max = iv[6], ci_min = iv[7], Topt_CASA = iv[8], sr_max = iv[9], sr_min = iv[10];
    // Loop through each n-days interval, accumulate the GPP
    var numDays = end_date.difference(start_date, 'day');
    var output = ee.List.sequence(0, numDays, interval)
        .map(function (i) {
            var start = start_date.advance(i, 'day');
            var end = start.advance(interval, 'day');
            var sv = calculate_shared_variables(start, end, co2_mm_gl);
            var T = sv[0], Tk = sv[1], Tdew = sv[2], Pmst = sv[3], sm = sv[4], actual_eva = sv[5], potential_eva = sv[6], vpd = sv[7], co2 = sv[8], lulc = sv[9], c3c4_map = sv[10], ci = sv[11];
            var par = calculate_par(par_param, start, end, interval);
            var fpar = calculate_fpar(fpar_param, start, end, sr_max, sr_min, interval);
            var lue_ci = calculate_lue_ci(lue_ci_param, c3c4_map, lulc, ci, ci_max, ci_min, resolution).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
            var ft = calculate_ft(ft_param, T, Tk, lulc, Topt_CASA, c3c4_map, resolution);
            var fw = calculate_fw(fw_param, start, end, T, Tdew, Pmst, sm, sm_max, sm_min, lswi_max, actual_eva, potential_eva, resolution);
            var fvpd_co2 = calculate_fvpd_co2(fvpd_co2_param, T, Tk, vpd, vpd_min, vpd_max, Pmst, co2, resolution, lulc);
            if (fw_param == '3-PG' && fvpd_co2_param == '3-PG') {
                var fw = fw.min(fvpd_co2);
                var fvpd_co2 = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
            } else if (ft_param == 'GLO-PEM | VPM | MVPM | EC-LUE-V3 | TEC | CCW' && fvpd_co2_param == 'EC-LUE-V3') {
                var ft = ee.Image(1.0).clipToBoundsAndScale({ 'geometry': geometry, 'scale': resolution });
            }
            // Calculate GPP
            var gpp = par.multiply(fpar).multiply(lue_ci).multiply(ft).multiply(fw).multiply(fvpd_co2);
            return par.addBands(fpar).addBands(lue_ci).addBands(ft).addBands(fw).addBands(fvpd_co2).addBands(gpp).select([0, 1, 2, 3, 4, 5, 6], ['PAR', 'fPAR', 'LUE_CI', 'fT', 'fW', 'fVPD_CO2', 'GPP']);
        });
    return ee.ImageCollection(output);
}

/*---------------------------------------------*/
// GUI code below
var par_param_list = ['MODIS18C2', 'BESS', 'ERA5-Land', 'MERRA-2', 'GLDAS-2.2', 'TerraClimate'];
var fpar_param_list = ['MOD15A2H', 'CASA', 'NDPI', 'GLO-PEM', 'C-Fix', 'VPM', 'EC-LUE-V3', 'TPM', 'Wang-model', 'LAI-based'];
var lue_ci_param_list = ['MOD17', 'P-model', 'C-Fix', 'CFLUX | CI-EF', 'Wang-model', 'VPM', 'TEC',  'MuSyQ', 'CCW', 'None'];
var ft_param_list = ['CASA', 'P-model', 'MOD17 | CFLUX | TL-LUE | CI-EF', 'C-Fix', 'GLO-PEM | VPM | MVPM | EC-LUE-V3 | TEC | CCW', 'CI-LUE', 'MuSyQ', 'Wang-model', 'None'];
var fw_param_list = ['P-model', 'GLO-PEM', '3-PG', 'VPM', 'TEC', 'CASA | MuSyQ', 'TCF', 'Wang-model', 'CI-EF', 'None'];
var fvpd_co2_param_list = ['P-model', '3-PG', 'MOD17 | CFLUX | TL-LUE', 'C-Fix', 'EC-LUE-V3', 'CCW', 'None'];
var resolution = 5000;
var geometry = ee.Geometry.Polygon({ coords: [[[-180, -90], [-180, 90], [180, 90], [180, -90]]], proj: "EPSG:4326", geodesic: false }); // global

// Create the UI elements.
var start_date_select = ui.DateSlider({ start: '2003-01-01', value: '2005-01-01' });
var end_date_select = ui.DateSlider({ start: '2003-01-01', value: '2006-01-01' });
var interval_select = ui.Slider({ min: 1, max: 30, value:10, step: 1 });
var resolution_select = ui.Slider({ min: 500, max: 5000, value:10, step: 100 });
var par_param_select = ui.Select({ items: par_param_list, value: par_param_list[0] });
var fpar_param_select = ui.Select({ items: fpar_param_list, value: fpar_param_list[0] });
var lue_ci_param_select = ui.Select({ items: lue_ci_param_list, value: lue_ci_param_list[0] });
var ft_param_select = ui.Select({ items: ft_param_list, value: ft_param_list[0] });
var fw_param_select = ui.Select({ items: fw_param_list, value: fw_param_list[0] });
var fvpd_co2_param_select = ui.Select({ items: fvpd_co2_param_list, value: fvpd_co2_param_list[0] });

// Add the UI elements to the map.
var parameter_panel = ui.Panel({
    widgets: [
        ui.Label('🙉 Start Date'),
        start_date_select,
        ui.Label('🐼 End Date'),
        end_date_select,
        // ui.Label('🐩 Resolution (m):'),
        // resolution_select,
        ui.Label('🦄 Interval Period (day):'),
        interval_select,
        ui.Label('🦬 Select PAR parameter:'),
        par_param_select,
        ui.Label('🐖 Select FPAR parameter:'),
        fpar_param_select,
        ui.Label('🐘 Select f(LUE,CI) parameter:'),
        lue_ci_param_select,
        ui.Label('🐁 Select f(T) parameter:'),
        ft_param_select,
        ui.Label('🐓 Select f(W) parameter:'),
        fw_param_select,
        ui.Label('🐍 Select f(VPD,CO2) parameter:'),
        fvpd_co2_param_select
    ],
    layout: 'flow',
    style: { fontWeight: 'bold', width: '250px', height: '750px' }
});

var inspectorPanel = ui.Panel({ style: { width: '30%'} });
var intro = ui.Panel([
    ui.Label({
        value: 'LUE fucntion values explorer',
        style: { fontSize: '20px', fontWeight: 'bold' }
    }),
    ui.Label('Click a location to see its time series of LUE fucntion values.')
]);
inspectorPanel.add(intro);
var lon = ui.Label();
var lat = ui.Label();
inspectorPanel.add(ui.Panel([lon, lat], ui.Panel.Layout.flow('horizontal')));

var generateChart = function (coords) {
    // Update the lon/lat panel with values from the click event.
    lon.setValue('lon: ' + coords.lon.toFixed(2));
    lat.setValue('lat: ' + coords.lat.toFixed(2));
    var point = ee.Geometry.Point(coords.lon, coords.lat);
    var dot = ui.Map.Layer(point, { color: '000000' }, 'clicked location');
    // Add the dot as the second layer, so it shows up on top of the composite.
    mapPanel.layers().set(1, dot);
    var result = main(ee.Date(start_date_select.getValue()[0]), ee.Date(end_date_select.getValue()[0]), resolution, interval_select.getValue(), par_param_select.getValue(), fpar_param_select.getValue(), lue_ci_param_select.getValue(), ft_param_select.getValue(), fw_param_select.getValue(), fvpd_co2_param_select.getValue());
    // Map over the collection to add the 'Index' property
    var resultWithId = result.map(function(img) {
        var index = ee.Number.parse(img.get('system:index')); // Convert 'system:index' to a number
        return img.set('Index', index); // Set the new 'time' property
    });
    // Sort the collection by the 'Index' property
    var result = resultWithId.sort('Index');
    var lue_param_chart = ui.Chart.image.series(result.select(['LUE_CI', 'fT', 'fW', 'fVPD_CO2', 'fPAR']), point, ee.Reducer.mean(), resolution, 'Index');
    var gpp_chart = ui.Chart.image.series(result.select(['PAR', 'GPP']), point, ee.Reducer.mean(), resolution, 'Index');
    lue_param_chart.setOptions({
        title: 'LUE fucntion values: time series',
        legend: { position: 'right' },
    });
    gpp_chart.setOptions({
        title: 'PAR(MJ d-1 m-2) & GPP(gC d-1  m-2) values: time series',
        legend: { position: 'right' },
    });
    // Add the chart at a fixed position, so that new charts overwrite older ones.
    inspectorPanel.widgets().set(2, lue_param_chart);
    inspectorPanel.widgets().set(3, gpp_chart);
};

var mapPanel = ui.Map();
mapPanel.onClick(generateChart);
mapPanel.style().set('cursor', 'crosshair');
// Configure the map.
ui.root.clear();
ui.root.add(parameter_panel);
ui.root.add(mapPanel);
ui.root.add(inspectorPanel);

"""
    SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes ice cover, surface/atmospheric/combined albedo, and net shortwave
flux from `Ts` and the current cloud climatology. Returns
`(SW, albedo, ice_cover)`.
"""
function SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Reuse workspace buffers
    ice_cover = ws.ice_cover_buf # output: ice fraction
    a_surf = ws.a_surf_buf       # surface albedo
    albedo = ws.albedo_buf       # output: combined albedo (surface + atmosphere)
    a_atmos = ws.a_atmos_buf     # atmospheric albedo
    sw = ws.sw_buf               # output: net shortwave flux

    z_topo = fields.z_topo
    glacier = fields.glacier

    # 1. Ice cover fraction
    @turbo for i in 1:xdim, j in 1:ydim
        T = Ts[i, j]
        land_expr = ifelse(T <= Tl_ice1, 1.0f0,
            ifelse(T < Tl_ice2,
                1.0f0 - (T - Tl_ice1) * inv_Tl_ice_range,
                0.0f0))
        ocean_expr = ifelse(T <= To_ice1, 1.0f0,
            ifelse(T < To_ice2,
                1.0f0 - (T - To_ice1) * inv_To_ice_range,
                0.0f0))
        ice_cover[i, j] = ifelse(z_topo[i, j] >= 0.0f0, land_expr, ocean_expr)
    end

    # 2. Atmospheric albedo
    cld = fields.cldclim
    ityr = timestate.ityr

    # 3. Surface albedo
    if cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            T = Ts[i, j]
            # Land albedo expression
            land_alb = ifelse(T <= Tl_ice1, a_no_ice + da_ice,
                ifelse(T >= Tl_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0f0 - (T - Tl_ice1) * inv_Tl_ice_range)))
            # Ocean albedo expression
            ocean_alb = ifelse(T <= To_ice1, a_no_ice + da_ice,
                ifelse(T >= To_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0f0 - (T - To_ice1) * inv_To_ice_range)))
            # Choose based on topography
            a_surf[i, j] = ifelse(z_topo[i, j] >= 0.0f0, land_alb, ocean_alb)
            # Glacier override: if glacier mask > 0.5, set to ice albedo
            a_surf[i, j] = ifelse(glacier[i, j] > 0.5f0, a_no_ice + da_ice, a_surf[i, j])
        end
    else
        @. a_surf = a_no_ice
    end

    # 4. albedo + shortwave flux.
    sw_solar = fields.sw_solar
    multiplier = state.sw_solar_forcing * 0.01f0 * S0_var
    @turbo for j in 1:ydim
        sf = sw_solar[j, ityr] * multiplier
        for i in 1:xdim
            aa = cld[i, j, ityr] * a_cloud
            a_atmos[i, j] = aa
            alb = a_surf[i, j] + aa - a_surf[i, j] * aa
            albedo[i, j] = alb
            sw[i, j] = sf * (1.0f0 - alb)
        end
    end

    return (SW=sw, albedo=albedo, ice_cover=ice_cover)
end

"""
    LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes atmospheric emissivity from CO₂/water-vapor/cloud columns, then
surface/upward/downward longwave flux. If `cfg.log_atmos_dmc` is false, only
`LW_down` is zeroed - `LW_up` is snapshotted beforehand and keeps its full
value (decouples surface from atmospheric downwelling feedback without
touching the atmosphere's own emission term). Returns
`(LW_surf, LW_up, LW_down, em)`.
"""
function LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Extract workspace buffers
    e_co2 = ws.e_co2_buf      # CO₂ [ppm scaled by pressure]
    e_vapor = ws.e_vapor_buf  # water vapour [kg/m²]
    em = ws.em_buf            # emissivity ε_atmos
    LW_surf = ws.LW_surf_buf  # surface long-wave flux [W/m²]
    LW_down = ws.LW_down_buf  # downward long-wave flux [W/m²]
    LW_up = ws.LW_up_buf      # upward long-wave flux [W/m²]

    wz_air = fields.wz_air
    co2_part = fields.co2_part
    ityr = timestate.ityr
    cldclim = fields.cldclim
    dTrad = fields.dTrad
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = p_emi

    # ── Effective columns, emissivity (log-regression, 10 parameters),
    # cloud adjustment, and surface/downward longwave flux
    @turbo for j in 1:ydim
        for i in 1:xdim
            e_vapor[i, j] = wz_air[i, j] * r_qviwv * q[i, j]
            e_co2[i, j] = wz_air[i, j] * CO2 * co2_part[i, j]
            em_val = p4 * log(p1 * e_co2[i, j] + p2 * e_vapor[i, j] + p3) +
                     p7 +
                     p5 * log(p1 * e_co2[i, j] + p3) +
                     p6 * log(p2 * e_vapor[i, j] + p3)
            em_val = (p8 - cldclim[i, j, ityr]) / p9 * (em_val - p10) + p10
            em[i, j] = em_val
            LW_surf[i, j] = -σ * Ts[i, j]^4
            LW_down_val = -em_val * σ * (Ta[i, j] + dTrad[i, j, ityr])^4
            LW_down[i, j] = LW_down_val
            LW_up[i, j] = LW_down_val
        end
    end

    if !cfg.log_atmos_dmc
        LW_down .= 0.0f0
    end

    return (LW_surf=LW_surf, LW_up=LW_up, LW_down=LW_down, em=em)
end

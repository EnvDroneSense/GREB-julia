# ── notebook cell 1df2b91b-be14-427a-87b3-95cdef26ce00  (orig lines 1359-1419) ──
"""
    SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes ice cover, surface/atmospheric/combined albedo, and net shortwave
flux from `Ts` and the current cloud climatology. Returns
`(SW, albedo, ice_cover)`.
"""
function SWradiation!(Ts, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Reuse workspace buffers
    ice_cover = ws.ice_cover_buf # output: ice fraction
    a_surf = ws.a_surf_buf    # surface albedo
    albedo = ws.albedo_buf    # output: combined albedo (surface + atmosphere) [hernoemen]
    a_atmos = ws.a_atmos_buf   # atmospheric albedo
    sw = ws.sw_buf        # output: net shortwave flux

    z_topo = fields.z_topo
    glacier = fields.glacier

    # 1. Ice cover fraction – branch‑free with ifelse, vectorized
    @turbo for i in 1:xdim, j in 1:ydim
        T = Ts[i, j]
        land_expr = ifelse(T <= Tl_ice1, 1.0,
            ifelse(T < Tl_ice2,
                1.0 - (T - Tl_ice1) * inv_Tl_ice_range,
                0.0))
        ocean_expr = ifelse(T <= To_ice1, 1.0,
            ifelse(T < To_ice2,
                1.0 - (T - To_ice1) * inv_To_ice_range,
                0.0))
        ice_cover[i, j] = ifelse(z_topo[i, j] >= 0.0, land_expr, ocean_expr)
    end

    # 2. Atmospheric albedo – simple multiplication, use broadcasting
    cld = @view fields.cldclim[:, :, timestate.ityr]
    @. a_atmos = cld * a_cloud

    # 3. Surface albedo – conditional logic, @turbo beneficial
    if cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            T = Ts[i, j]
            # Land albedo expression
            land_alb = ifelse(T <= Tl_ice1, a_no_ice + da_ice,
                ifelse(T >= Tl_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0 - (T - Tl_ice1) * inv_Tl_ice_range)))
            # Ocean albedo expression
            ocean_alb = ifelse(T <= To_ice1, a_no_ice + da_ice,
                ifelse(T >= To_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0 - (T - To_ice1) * inv_To_ice_range)))
            # Choose based on topography
            a_surf[i, j] = ifelse(z_topo[i, j] >= 0.0, land_alb, ocean_alb)
            # Glacier override: if glacier mask > 0.5, set to ice albedo
            a_surf[i, j] = ifelse(glacier[i, j] > 0.5, a_no_ice + da_ice, a_surf[i, j])
        end
    else
        @. a_surf = a_no_ice
    end

    # 4. Combined albedo
    @. albedo = a_surf + a_atmos - a_surf * a_atmos

    # 5. Shortwave flux
    multiplier = state.sw_solar_forcing * 0.01 * S0_var
    sw_solar = fields.sw_solar
    for j in 1:ydim
        sf = sw_solar[j, timestate.ityr] * multiplier
        # @views: `albedo[:, j]` on the RHS of `@.` is plain getindex (not a
        # dotview like the LHS), so without this it allocates a fresh
        # Vector{Float64} every iteration.
        @views @. sw[:, j] = sf * (1.0 - albedo[:, j])
    end

    return (SW=sw, albedo=albedo, ice_cover=ice_cover)
end

# ── notebook cell c0c40037-4169-4d38-bebe-2086cebc24f2  (orig lines 1437-1476) ──
"""
    LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes atmospheric emissivity from CO₂/water-vapor/cloud columns, then
surface/upward/downward longwave flux. If `cfg.log_atmos_dmc` is false, only
`LW_down` is zeroed — `LW_up` is snapshotted beforehand and keeps its full
value (decouples surface from atmospheric downwelling feedback without
touching the atmosphere's own emission term). Returns
`(LW_surf, LW_up, LW_down, em)`.
"""
function LWradiation!(Ts, Ta, q, CO2, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Extract workspace buffers
    e_co2 = ws.e_co2_buf    # CO₂ [ppm scaled by pressure]
    e_vapor = ws.e_vapor_buf  # water vapour [kg/m²]
    em = ws.em_buf       # emissivity ε_atmos
    LW_surf = ws.LW_surf_buf  # surface long-wave flux [W/m²]
    LW_down = ws.LW_down_buf  # downward long-wave flux [W/m²]
    LW_up = ws.LW_up_buf    # upward long-wave flux [W/m²]

    wz_air = fields.wz_air
    co2_part = fields.co2_part

    # Current cloud cover (climatology, 3D array)
    e_cloud = @view fields.cldclim[:, :, timestate.ityr]

    ## 1. Effective columns (topography scaling via wz_air, precomputed)
    @. e_vapor = wz_air * r_qviwv * q
    @. e_co2 = wz_air * CO2 * co2_part

    # ── Emissivity (log-regression with 10 parameters) ─────────
    @. em = p_emi[4] * log(p_emi[1] * e_co2 + p_emi[2] * e_vapor + p_emi[3]) +
            p_emi[7] +
            p_emi[5] * log(p_emi[1] * e_co2 + p_emi[3]) +
            p_emi[6] * log(p_emi[2] * e_vapor + p_emi[3])

    # Cloud adjustment
    @. em = (p_emi[8] - e_cloud) / p_emi[9] * (em - p_emi[10]) + p_emi[10]

    # 4. Radiation temperature (precomputed offset dTrad = -0.16*Tclim - 5 K)
    dTr = @view fields.dTrad[:, :, timestate.ityr]
    @. LW_surf = -σ * Ts^4
    @. LW_down = -em * σ * (Ta + dTr)^4
    LW_up .= LW_down

    if !cfg.log_atmos_dmc
        LW_down .= 0.0
    end

    return (LW_surf=LW_surf, LW_up=LW_up, LW_down=LW_down, em=em)
end

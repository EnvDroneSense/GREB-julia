"""
    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)

Runs one timestep's physics pipeline in order — [`SWradiation!`](@ref) →
[`LWradiation!`](@ref) → sensible heat → [`hydro!`](@ref) →
[`circulation!`](@ref) (temperature, then humidity) → [`deep_ocean!`](@ref) —
and returns a named tuple of every intermediate flux/tendency needed by
[`diagnostics!`](@ref) and the caller's own state update.
"""
# ── notebook cell e493fae7-239a-494c-9a59-728446d70f7a  (orig lines 2080-2127) ──
function tendencies!(CO2, Ts, Ta, To, q, fields::ClimateFields, state::ModelState, ws::CirculationWorkspace,
    timestate, cfg::PhysicsConfig)

    # Short-wave radiation → albedo, SW flux
    sw_out = SWradiation!(Ts, fields, state, timestate, cfg, ws)

    # Long-wave radiation → LW_surf, LW_up, LW_down, emissivity
    lw_out = LWradiation!(Ts, Ta, q, CO2, fields, timestate, cfg, ws)

    # Sensible heat flux
    Q_sens = ws.Q_sens_buf
    if cfg.log_atmos_dmc
        @. Q_sens = ct_sens * (Ta - Ts)
    else
        fill!(Q_sens, 0.0)
    end

    # Hydrological cycle → latent heat + evaporation/rain tendencies
    hy_out = hydro!(Ts, q, fields, timestate, cfg, ws)

    # Atmospheric circulation — temperature diffusion/advection
    circulation!(Ta, z_air, ws.dTa_crcl, fields, ws, timestate, cfg)

    # Atmospheric circulation — water-vapour diffusion/advection
    circulation!(q, z_vapor, ws.dq_crcl, fields, ws, timestate, cfg)

    # Deep ocean coupling
    do_out = deep_ocean!(Ts, To, fields, timestate, cfg, ws)

    return (albedo=sw_out.albedo,
        SW=sw_out.SW,
        ice_cover=sw_out.ice_cover,
        LW_surf=lw_out.LW_surf,
        Q_lat=hy_out.Q_lat,
        Q_sens=Q_sens,
        Q_lat_air=hy_out.Q_lat_air,
        dq_eva=hy_out.dq_eva,
        dq_rain=hy_out.dq_rain,
        dq_crcl=ws.dq_crcl,
        dTa_crcl=ws.dTa_crcl,
        dT_ocean=do_out.dT_ocean,
        dTo=do_out.dTo,
        LW_down=lw_out.LW_down,
        LW_up=lw_out.LW_up,
        em=lw_out.em)
end

"""
    forcing(it, year, cfg::PhysicsConfig, fields::ClimateFields, icmn_ctrl; nstep_yr=nstep_yr)

Returns `(CO2, sw_solar_forcing)` for the current timestep, computed
according to `cfg.experiment`. Some experiments (the `regional_co2_*` family)
also mutate `fields.co2_part` as a side effect. `:full_model` short-circuits
before the experiment dispatch chain.
"""
# ── notebook cell 1894ad94-cdf8-4e79-a0e5-b72088db31be  (orig lines 2162-2343) ──
function forcing(it, year, cfg::PhysicsConfig, fields::ClimateFields, icmn_ctrl; nstep_yr=nstep_yr)
    # Default CO₂ concentration
    CO2 = cfg.co2_concentration
    sw_solar_forcing = 1.0

    # Fast path for the main experiment: every branch below is a no-op for
    # :full_model, but without this it still falls through all of them,
    # including an allocating `string(cfg.experiment)` + `startswith` check
    # (line ~154) on every single timestep of the most common run type.
    if cfg.experiment == :full_model
        return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
    end

    # 📜 Legacy experiments ───────────
    if cfg.experiment == :constant_topo
        CO2 = 550.0  # 550 ppm CO₂ steady state

    elseif cfg.experiment == :a1b_scenario
        CO2_1950 = 310.0;
        CO2_2000 = 370.0;
        CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end

        # 💨 CO₂ scaling experiments ──────────────────────────────────────────────
    elseif cfg.experiment == :co2_double
        CO2 = 680.0  # 2×CO₂ (already set, but explicit)

    elseif cfg.experiment == :co2_quadruple
        CO2 = 1360.0  # 4×CO₂

    elseif cfg.experiment == :co2_10x
        CO2 = 3400.0  # 10×CO₂

    elseif cfg.experiment == :co2_half
        CO2 = 170.0  # 0.5×CO₂

    elseif cfg.experiment == :co2_zero
        CO2 = 0.0  # 0×CO₂ (no greenhouse effect)

    # ☀️ Solar forcing experiments ───────────────────────────────────────────
    elseif cfg.experiment == :solar_plus27
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 27.0) / 1365.0

    elseif cfg.experiment == :solar_cycle_11yr
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 1.0 * sin(2π * year / 11.0)) / 1365.0

        # 📈 Enhanced A1B scenario ──────────────────────────────────────────────
    elseif cfg.experiment == :a1b_enhanced
        CO2_1950 = 310.0;
        CO2_2000 = 370.0;
        CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end

        # ── Time-varying CO₂ experiments ────────────
    elseif cfg.experiment == :co2_sine_wave
        CO2 = 340.0 + 170.0 + 170.0 * cos(2π * (year - 13.0) / 30.0)

    elseif cfg.experiment == :co2_step
        CO2 = year >= 1980 ? 340.0 : 680.0

        # ── Paleoclimate experiments ────────────────────
    elseif cfg.experiment == :paleo_231kyr
        CO2 = 200.0

    elseif cfg.experiment == :paleo_solar_modern_co2
        CO2 = 340.0

    elseif cfg.experiment == :modern_solar_paleo_co2
        CO2 = 200.0

        # ── Orbital forcing experiments ─────────────────
    elseif cfg.experiment == :obliquity
        CO2 = 340.0     # Solar forcing loaded externally

    elseif cfg.experiment == :eccentricity
        CO2 = 340.0     # Solar forcing loaded externally

    elseif cfg.experiment == :earth_sun_distance
        CO2 = 340.0     # Solar constant varies with Earth-Sun distance
        sw_solar_forcing = (1.0 / (1.0 + 0.01 * cfg.earth_sun_distance_pct))^2

    # 📂 File I/O dependent experiments (placeholders) ───────────────────────
    elseif cfg.experiment == :rcp26
        error("RCP2.6 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp45
        error("RCP4.5 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp60
        error("RCP6.0 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp85
        CO2 = 340.0  # Handled by boundary conditions

    elseif cfg.experiment == :custom_co2
        error("Custom CO₂ scenario requires external trajectory file. Not yet implemented.")

        # 🌍 Regional/partial CO₂ experiments — static masks ─────────────────────
        # nh/sh/tropics/extratropics depend only on fixed latitude-row ranges,
        # so their mask is set once in init_model! (IMPROVEMENTS.md §4.4) —
        # nothing to recompute here every timestep.
    elseif cfg.experiment in (:regional_co2_nh, :regional_co2_sh, :regional_co2_tropics, :regional_co2_extratropics)
        CO2 = 680.0

        # 🌍 Regional/partial CO₂ experiments — dynamic masks ────────────────────
        # ocean/land_ice depend on icmn_ctrl (the control run's ice climatology),
        # which doesn't exist at init_model! time — these two must stay
        # recomputed every timestep. winter/summer don't touch co2_part at all
        # (season-dependent, not region-dependent).
    elseif startswith(string(cfg.experiment), "regional_co2_")
        co2_part = fields.co2_part
        co2_part .= 1.0

        if cfg.experiment == :regional_co2_ocean
            # 2×CO₂ Ocean only
            CO2 = 680.0
            z_topo = fields.z_topo
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] > 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 0.5
                end
            end

        elseif cfg.experiment == :regional_co2_land_ice
            # 2×CO₂ Land/Ice only
            CO2 = 680.0
            z_topo = fields.z_topo
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] <= 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 1.0
                end
            end

        elseif cfg.experiment == :regional_co2_winter
            # 2×CO₂ Boreal Winter only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 680.0 : 340.0

        elseif cfg.experiment == :regional_co2_summer
            # 2×CO₂ Boreal Summer only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 340.0 : 680.0
        end

        # 📊 Forced boundary condition experiments (handled in scenario loop) ─────
    elseif cfg.experiment == :elnino || cfg.experiment == :lanina || cfg.experiment == :rcp85
        CO2 = 340.0
    end

    return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
end

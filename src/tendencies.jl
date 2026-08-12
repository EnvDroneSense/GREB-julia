"""
    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg; ws_a=ws, ws_q=ws)

Runs one timestep's physics pipeline - [`SWradiation!`](@ref) →
[`LWradiation!`](@ref) → sensible heat → [`hydro!`](@ref) →
[`circulation!`](@ref) (temperature, then humidity) → [`deep_ocean!`](@ref) -
and returns a named tuple of every intermediate flux/tendency needed by
[`diagnostics!`](@ref) and the caller's own state update.

The two `circulation!` calls are independent of each other and of every
other stage (each reads only pre-timestep state and writes disjoint
buffers), so when the caller supplies distinct `ws_a`/`ws_q` workspaces
*and* `Threads.nthreads() > 1`, they run concurrently via `Threads.@spawn`
while the remaining stages run on `ws`. With the default `ws_a=ws_q=ws`
"""
function tendencies!(CO2, Ts, Ta, To, q, fields::ClimateFields, state::ModelState, ws::CirculationWorkspace,
    timestate, cfg::PhysicsConfig; ws_a::CirculationWorkspace=ws, ws_q::CirculationWorkspace=ws)

    parallel = Threads.nthreads() > 1 && ws_a !== ws_q

    # Atmospheric circulation - temperature/water-vapour diffusion/advection.
    if parallel
        t_a = Threads.@spawn circulation!(Ta, z_air, ws_a.dTa_crcl, fields, ws_a, timestate, cfg)
        t_q = Threads.@spawn circulation!(q, z_vapor, ws_q.dq_crcl, fields, ws_q, timestate, cfg)
    else
        circulation!(Ta, z_air, ws_a.dTa_crcl, fields, ws_a, timestate, cfg)
        circulation!(q, z_vapor, ws_q.dq_crcl, fields, ws_q, timestate, cfg)
    end

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

    # Deep ocean coupling
    do_out = deep_ocean!(Ts, To, fields, timestate, cfg, ws)

    if parallel
        wait(t_a)
        wait(t_q)
    end

    return (albedo=sw_out.albedo,
        SW=sw_out.SW,
        ice_cover=sw_out.ice_cover,
        LW_surf=lw_out.LW_surf,
        Q_lat=hy_out.Q_lat,
        Q_sens=Q_sens,
        Q_lat_air=hy_out.Q_lat_air,
        dq_eva=hy_out.dq_eva,
        dq_rain=hy_out.dq_rain,
        dq_crcl=ws_q.dq_crcl,
        dTa_crcl=ws_a.dTa_crcl,
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
before the experiment dispatch chain. The `:ssp*`/`:historical_co2`
experiments look `year` up in `cfg.co2_scenario`.
"""
function forcing(it, year, cfg::PhysicsConfig, fields::ClimateFields, icmn_ctrl; nstep_yr=nstep_yr)
    # Default CO₂ concentration
    CO2 = cfg.co2_concentration
    sw_solar_forcing = 1.0

    # Fast path for the main experiment
    if cfg.experiment == :full_model
        return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
    end

    # - Legacy experiments ───────────
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

    # - CO₂ scaling experiments ──────────────────────────────────────────────
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

    # - Solar forcing experiments ───────────────────────────────────────────
    elseif cfg.experiment == :solar_plus27
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 27.0) / 1365.0

    elseif cfg.experiment == :solar_cycle_11yr
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 1.0 * sin(2π * year / 11.0)) / 1365.0

    # - Enhanced A1B scenario ──────────────────────────────────────────────
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

    # - File I/O dependent experiments (placeholders) ───────────────────────
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

    # - IPCC SSP/historical scenarios - CO₂ read from a per-year lookup table ─
    elseif cfg.experiment in (:ssp119, :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2)
        yr = round(Int, year)
        haskey(cfg.co2_scenario, yr) ||
            error("No CO2 data for year $yr in $(cfg.experiment) scenario table " *
                  "(loaded $(length(cfg.co2_scenario)) years)")
        CO2 = cfg.co2_scenario[yr]

    # - Regional/partial CO₂ experiments - static masks ─────────────────────
    elseif cfg.experiment in (:regional_co2_nh, :regional_co2_sh, :regional_co2_tropics, :regional_co2_extratropics)
        CO2 = 680.0

    # - Regional/partial CO₂ experiments - dynamic masks ────────────────────
    # `:regional_co2_ocean`/`:regional_co2_land_ice`'s mask depends only on
    # `z_topo` and `icmn_ctrl` (both fixed for the whole scenario run), so
    # it's computed once (`it == 1`) rather than recomputed every timestep;
    # `fields.co2_part` holds the correct static mask for every later
    # timestep since nothing else touches it for these two experiments.
    elseif startswith(string(cfg.experiment), "regional_co2_")
        if cfg.experiment == :regional_co2_ocean
            # 2×CO₂ Ocean only
            CO2 = 680.0
            if it == 1
                co2_part = fields.co2_part
                co2_part .= 1.0
                z_topo = fields.z_topo
                for j in 1:ydim, i in 1:xdim
                    if z_topo[i, j] > 0.0
                        co2_part[i, j] = 0.5
                    end
                end
                # Annual-mean ice cover
                icmn_ctrl1 = dropdims(sum(icmn_ctrl, dims=3), dims=3) ./ size(icmn_ctrl, 3)
                for j in 1:ydim, i in 1:xdim
                    if icmn_ctrl1[i, j] >= 0.5
                        co2_part[i, j] = 0.5
                    end
                end
            end

        elseif cfg.experiment == :regional_co2_land_ice
            # 2×CO₂ Land/Ice only
            CO2 = 680.0
            if it == 1
                co2_part = fields.co2_part
                co2_part .= 1.0
                z_topo = fields.z_topo
                for j in 1:ydim, i in 1:xdim
                    if z_topo[i, j] <= 0.0
                        co2_part[i, j] = 0.5
                    end
                end
                icmn_ctrl1 = dropdims(sum(icmn_ctrl, dims=3), dims=3) ./ size(icmn_ctrl, 3)
                for j in 1:ydim, i in 1:xdim
                    if icmn_ctrl1[i, j] >= 0.5
                        co2_part[i, j] = 1.0
                    end
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

        # - Forced boundary condition experiments (handled in scenario loop) ─────
    elseif cfg.experiment == :elnino || cfg.experiment == :lanina || cfg.experiment == :rcp85
        CO2 = 340.0
    end

    return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
end

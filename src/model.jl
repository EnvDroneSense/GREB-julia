# ── notebook cell d404043f-8080-4262-9ab6-d9bb13eee504  (orig lines 1200-1304) ──
function init_model!(cfg::PhysicsConfig)

    # ── Hydrology Parameter Initialization ────────────
    set_hydrology_parameters!(cfg)

    # ── Reset regional CO₂ mask ───────────────────────
    # co2_part is a module global only ever mutated by forcing()'s
    # regional_co2_* branches; without this reset a regional experiment
    # run earlier in the same session would leak its stale mask into any
    # later, unrelated run.
    co2_part .= 1.0

    # ── dTrad: offset between T_atm and radiation temperature ────
    @. dTrad = -0.16 * Tclim - 5.0

    # ── z_ocean: 3× maximum mixed-layer depth over the year ──────
    z_ocean .= 3.0 .* dropdims(maximum(mldclim; dims=3); dims=3)

    # ── Sensitivity experiment overrides ─────────────────────────
    if !cfg.log_topo_drsp
        @. z_topo = min(z_topo, 1.0)       # constant topography
    end

    # Apply cloud deconstruction switch
    if !cfg.log_clouds_dmc
        cldclim .= 0.0  # zero cloud climatology
    end

    # Apply flux correction conditional zeroing (MSCM feature)
    if !cfg.log_topo_drsp && !cfg.log_qflux_dmc
        TF_correct .= 0.0
        qF_correct .= 0.0
        ToF_correct .= 0.0
    end

    # Climatology modifications
    if !cfg.log_hydro_dmc
        qclim .= 0.0  # zero out humidity climatology
    end

    if !cfg.log_vapor_dmc
        qclim .= 0.0052          # constant water vapor
    end

    if !cfg.log_ocean_dmc
        mldclim .= d_ocean       # no deep ocean
    end

    # ── Experiment Handler ───────────────────────────────────
    # Apply advanced experiment forcing
    if cfg.experiment == :rcp85
        @info "Applying CMIP5 RCP8.5 climate change forcing"
        Tclim .+= Tclim_anom_cc
        uclim .+= uclim_anom_cc
        vclim .+= vclim_anom_cc
        omegaclim .+= omegaclim_anom_cc
        wsclim .+= wsclim_anom_cc
    elseif cfg.experiment == :elnino
        @info "Applying ERA-Interim El Niño forcing"
        Tclim .+= Tclim_anom_enso
        uclim .+= uclim_anom_enso
        vclim .+= vclim_anom_enso
        omegaclim .+= omegaclim_anom_enso
        wsclim .+= wsclim_anom_enso
    elseif cfg.experiment == :lanina
        @info "Applying ERA-Interim La Niña forcing"
        Tclim .-= Tclim_anom_enso
        uclim .-= uclim_anom_enso
        vclim .-= vclim_anom_enso
        omegaclim .-= omegaclim_anom_enso
        wsclim .-= wsclim_anom_enso
    end

    # ── Topography pressure weights ─────
    @. wz_air = exp(-z_topo / z_air)
    @. wz_vapor = exp(-z_topo / z_vapor)

    # ── Surface heat capacity ────────────────────────────────────
    @inbounds for j in 1:ydim
        for i in 1:xdim
            if z_topo[i, j] > 0.0
                cap_surf[i, j] = cap_land
            else
                cap_surf[i, j] = cfg.log_ocean_dmc ? cap_ocean * mldclim[i, j, 1] : cap_land
            end
        end
    end

    # ── Initial conditions from last time step of climatology ────
    Ts_ini = Tclim[:, :, nstep_yr] |> copy   # surface temperature
    Ta_ini = copy(Ts_ini)                      # air temperature = Tsurf
    To_ini = Toclim[:, :, nstep_yr] |> copy   # deep ocean temperature
    q_ini = qclim[:, :, nstep_yr] |> copy   # atmospheric water vapor

    # ── Control CO₂ level ───────────────────────────────────────
    CO2_ctrl = cfg.co2_concentration

    if cfg.experiment == :a1b_scenario
        CO2_ctrl = 298.0  # A1B scenario baseline
    elseif cfg.experiment in (:a1b_enhanced, :rcp26, :rcp45, :rcp60, :rcp85, :custom_co2)
        CO2_ctrl = 280.0  # IPCC scenarios baseline
    end

    if !cfg.log_co2_dmc
        CO2_ctrl = 0.0  # Zero CO2 for deconstruction experiments
    end

    return (Ts_ini=Ts_ini, Ta_ini=Ta_ini, To_ini=To_ini,
        q_ini=q_ini, CO2_ctrl=CO2_ctrl)
end

# ── notebook cell 584e767e-4dc5-4821-af63-d6a825326d9e  (orig lines 2869-2928) ──
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace, time_flux)
    for it in 1:(time_flux*ndt_days*ndays_yr)
        timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
        timestate.ityr = mod(it - 1, nstep_yr) + 1
        ityr = timestate.ityr

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, ws, timestate, cfg)

        # Views into climatology & correction fields
        Tc = @view Tclim[:, :, ityr]
        Toc = @view Toclim[:, :, ityr]
        qc = @view qclim[:, :, ityr]
        TFc = @view TF_correct[:, :, ityr]
        ToFc = @view ToF_correct[:, :, ityr]
        qFc = @view qF_correct[:, :, ityr]

        # ── Surface temperature ──────────────────────────────
        # Uncorrected state (store in workspace buffer)
        @. ws.Ts0_buf = Ts + tend.dT_ocean + Δt * (
            tend.SW + tend.LW_surf - tend.LW_down +
            tend.Q_lat + tend.Q_sens
        ) / cap_surf

        # Correction and corrected state
        @. TFc = (Tc - ws.Ts0_buf) * cap_surf / Δt
        @. ws.Ts0_buf = ws.Ts0_buf + TFc * Δt / cap_surf

        # ── Air temperature ──────────────────────────────────
        @. ws.Ta0_buf = Ta + tend.dTa_crcl + ΔT_AIR_FACTOR * (
            tend.LW_up + tend.LW_down - tend.em * tend.LW_surf +
            tend.Q_lat_air - tend.Q_sens
        )

        # ── Deep ocean ───────────────────────────────────────
        @. ws.To0_buf = To + tend.dTo
        @. ToFc = Toc - ws.To0_buf
        @. ws.To0_buf = ws.To0_buf + ToFc

        # ── Humidity ─────────────────────────────────────────
        @. ws.q0_buf = q + tend.dq_crcl + Δt * (tend.dq_eva + tend.dq_rain)
        @. qFc = qc - ws.q0_buf
        @. ws.q0_buf = ws.q0_buf + qFc

        # Sea ice (updates cap_surf in place)
        seaice!(ws.Ts0_buf, timestate, cfg)

        # Diagnostics
        diagnostics!(it, 0.0, CO2_ctrl, ws.Ts0_buf, ws.Ta0_buf, ws.To0_buf, ws.q0_buf,
            tend.albedo, tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

        # Advance state (broadcast – same as `.=`)
        @. Ts = ws.Ts0_buf
        @. Ta = ws.Ta0_buf
        @. q = ws.q0_buf
        @. To = ws.To0_buf
    end
    return nothing
end

# ── notebook cell 92c3bd68-bd07-4381-9c04-e6611650cd1e  (orig lines 2929-3043) ──
function greb_model!(time_flux, time_ctrl, time_scnr, cfg::PhysicsConfig; jld2_dir::AbstractString="")

    # ── 1. Initialisation ───────────────────────────────────────
    ini = init_model!(cfg)
    Ts_ini = ini.Ts_ini;
    Ta_ini = ini.Ta_ini
    To_ini = ini.To_ini;
    q_ini = ini.q_ini
    CO2_ctrl = ini.CO2_ctrl

    # Determine experiment type from cfg
    is_orbital_exp = cfg.experiment in (:obliquity, :eccentricity, :earth_sun_distance)
    is_forced_boundary = cfg.experiment in (:rcp85, :elnino, :lanina)
    is_sst_plus1 = cfg.experiment == :sst_plus1

    # Workspace and accumulator
    ws = CirculationWorkspace()
    acc = MonthlyAccumulator()

    # Initialize time state
    timestate = TimeState(1, 1)

    # ── 2. Flux-correction spin-up ──────────────────────────────
    if cfg.log_topo_drsp || cfg.log_qflux_dmc
        if !cfg.log_topo_drsp && cfg.log_qflux_dmc
            println("% loading flux correction fields...")
            load_flux_corrections_jld2!(jld2_dir)
        end
        println("% flux correction  CO2 = ", CO2_ctrl)
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, timestate, cfg, ws, time_flux)
    else
        println("Flux correction skipped")
    end

    # Reset accumulators after spin-up
    reset!(acc)

    # ── 3. Control run ──────────────────────────────────────────
    println("CONTROL RUN: CO2 = ", CO2_ctrl, " time = ", time_ctrl, " yr")

    # Initialize state arrays
    Ts = copy(Ts_ini);
    Ta = copy(Ta_ini)
    To = copy(To_ini);
    q = copy(q_ini)
    sw_solar_forcing_state[] = 1.0
    mon = 1;
    year = 1970;
    irec = 0

    ctrl_output = MonthlyRecord[]
    sizehint!(ctrl_output, time_ctrl * 12)  # Pre-allocate for all months
    timestate = TimeState(1, 1)  # Initialize time state

    for it in 1:(time_ctrl*nstep_yr)
        (mon, irec) = time_loop!(it, year, CO2_ctrl, mon, irec,
            Ts, Ta, q, To, ctrl_output, ws, acc, timestate, cfg)
        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

    # ── Build ice climatology from control output ───────────────
    # Compute annual mean ice cover from the stored control monthly means
    # (Assuming control output contains monthly ice cover; adjust if needed)
    ice_forcing = compute_annual_ice_climatology(ctrl_output)

    # ── 4. Scenario run ─────────────────────────────────────────
    println("SCENARIO: ", cfg.experiment, "  time = ", time_scnr, " yr")

    # Reset state to initial conditions
    Ts .= Ts_ini;
    Ta .= Ta_ini
    q .= q_ini;
    To .= To_ini
    year = is_orbital_exp ? 1 : 1950
    CO2 = 340.0;
    mon = 1;
    irec = 0

    sw_solar_forcing_state[] = 1.0
    reset!(acc)  # Use accumulator reset

    scnr_output = MonthlyRecord[]
    if time_scnr > 0
        sizehint!(scnr_output, time_scnr * 12)
    end

    for it in 1:(time_scnr*nstep_yr)
        # Obtain forcing (CO2 and solar multiplier)
        forcing_result = forcing(it, year, cfg, ice_forcing; nstep_yr=nstep_yr)
        CO2 = forcing_result.CO2
        sw_solar_forcing_state[] = forcing_result.sw_solar_forcing

        # Forced‑boundary experiments: overwrite Ts with climatology
        if is_forced_boundary
            ityr_now = mod(it - 1, nstep_yr) + 1
            Ts .= @view Tclim[:, :, ityr_now]
        end

        # SST+1 K experiment
        if is_sst_plus1
            CO2 = CO2_ctrl
            ityr_now = mod(it - 1, nstep_yr) + 1
            @. Ts = ifelse(z_topo < 0.0, Tclim[:, :, ityr_now] + 1.0, Ts)
        end

        (mon, irec) = time_loop!(it, year, CO2, mon, irec,
            Ts, Ta, q, To, scnr_output, ws, acc, timestate, cfg)

        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

    # Post‑processing: anomalies for non‑orbital experiments
    if !is_orbital_exp && !isempty(ctrl_output) && !isempty(scnr_output)
        ctrl_clim = build_monthly_climatology(ctrl_output)
        scnr_output = apply_scenario_anomalies(scnr_output, ctrl_clim)
    end

    return (ctrl=ctrl_output, scnr=scnr_output)
end

"""
    init_model!(cfg::PhysicsConfig, fields::ClimateFields)

One-time per-run setup: derives `cfg`'s hydrology parameters, resets the
regional-CO₂ mask, applies CO₂-response climatology overrides
(`log_clouds_drsp`/`log_humid_drsp`/`log_ocean_drsp`), and computes the
control-run initial state. Returns `(Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)`.
"""
function init_model!(cfg::PhysicsConfig, fields::ClimateFields)

    # ── Hydrology Parameter Initialization ────────────
    set_hydrology_parameters!(cfg)

    # ── Reset regional CO₂ mask ───────────────────────
    # Only matters if `fields` is being reused across multiple greb_model!
    # calls (the default is a fresh ClimateFields() per call); a regional
    # experiment run earlier against the same fields instance would otherwise
    # leak its stale mask into a later, unrelated run against it.
    fields.co2_part .= 1.0f0

    if cfg.experiment == :regional_co2_nh
        fields.co2_part[:, 1:24] .= 0.5f0
    elseif cfg.experiment == :regional_co2_sh
        fields.co2_part[:, 25:48] .= 0.5f0
    elseif cfg.experiment == :regional_co2_tropics
        fields.co2_part[:, 1:15] .= 0.5f0
        fields.co2_part[:, 33:48] .= 0.5f0
        for i in 4:4:96
            fields.co2_part[i, 33] = 1.0f0
            fields.co2_part[i, 15] = 1.0f0
        end
    elseif cfg.experiment == :regional_co2_extratropics
        fields.co2_part[:, 16:32] .= 0.5f0
        for i in 4:4:96
            fields.co2_part[i, 32] = 1.0f0
            fields.co2_part[i, 16] = 1.0f0
        end
    end

    Tclim = fields.Tclim
    z_topo = fields.z_topo

    # ── dTrad: offset between T_atm and radiation temperature ────
    @. fields.dTrad = -0.16f0 * Tclim - 5.0f0

    # ── z_ocean: 3× maximum mixed-layer depth over the year ──────
    fields.z_ocean .= 3.0f0 .* dropdims(maximum(fields.mldclim; dims=3); dims=3)

    # ── Sensitivity experiment overrides ─────────────────────────
    if !cfg.log_topo_drsp
        @. z_topo = min(z_topo, 1.0f0)       # constant topography
    end

    # Apply cloud deconstruction switch
    if !cfg.log_clouds_dmc
        fields.cldclim .= 0.0f0  # zero cloud climatology
    end

    # Apply flux correction conditional zeroing (MSCM feature)
    if !cfg.log_topo_drsp && !cfg.log_qflux_dmc
        fields.TF_correct .= 0.0f0
        fields.qF_correct .= 0.0f0
        fields.ToF_correct .= 0.0f0
    end

    # Climatology modifications
    if !cfg.log_hydro_dmc
        fields.qclim .= 0.0f0  # zero out humidity climatology
    end

    if !cfg.log_clouds_drsp
        fields.cldclim .= 0.7f0           # constant cloud cover (2xCO2 deconstruction)
    end

    if !cfg.log_humid_drsp
        fields.qclim .= 0.0052f0          # constant water vapor
    end

    if !cfg.log_ocean_drsp
        fields.mldclim .= d_ocean       # no deep ocean
    end

    # ── Experiment Handler ───────────────────────────────────
    # Apply advanced experiment forcing
    if cfg.experiment == :rcp85
        @info "Applying CMIP5 RCP8.5 climate change forcing"
        Tclim .+= fields.Tclim_anom_cc
        fields.uclim .+= fields.uclim_anom_cc
        fields.vclim .+= fields.vclim_anom_cc
        fields.omegaclim .+= fields.omegaclim_anom_cc
        fields.wsclim .+= fields.wsclim_anom_cc
    elseif cfg.experiment == :elnino
        @info "Applying ERA-Interim El Niño forcing"
        Tclim .+= fields.Tclim_anom_enso
        fields.uclim .+= fields.uclim_anom_enso
        fields.vclim .+= fields.vclim_anom_enso
        fields.omegaclim .+= fields.omegaclim_anom_enso
        fields.wsclim .+= fields.wsclim_anom_enso
    elseif cfg.experiment == :lanina
        @info "Applying ERA-Interim La Niña forcing"
        Tclim .-= fields.Tclim_anom_enso
        fields.uclim .-= fields.uclim_anom_enso
        fields.vclim .-= fields.vclim_anom_enso
        fields.omegaclim .-= fields.omegaclim_anom_enso
        fields.wsclim .-= fields.wsclim_anom_enso
    end

    # ── Topography pressure weights ─────
    @. fields.wz_air = exp(-z_topo / z_air)
    @. fields.wz_vapor = exp(-z_topo / z_vapor)

    # ── hydro!'s log_rain==1 rain-limit divisor: depends only on wz_vapor,
    # invariant for the whole run
    @. fields.rain_limit = -0.0015f0 / (fields.wz_vapor * r_qviwv * 86400.0f0)

    # ── Surface heat capacity ────────────────────────────────────
    cap_surf = fields.cap_surf
    mldclim = fields.mldclim
    @inbounds for j in 1:ydim
        for i in 1:xdim
            if z_topo[i, j] > 0.0f0
                cap_surf[i, j] = cap_land
            else
                cap_surf[i, j] = cfg.log_ocean_dmc ? cap_ocean * mldclim[i, j, 1] : cap_land
            end
        end
    end

    # ── Initial conditions from last time step of climatology ────
    Ts_ini = Tclim[:, :, nstep_yr] |> copy          # surface temperature
    Ta_ini = copy(Ts_ini)                           # air temperature = Tsurf
    To_ini = fields.Toclim[:, :, nstep_yr] |> copy  # deep ocean temperature
    q_ini = fields.qclim[:, :, nstep_yr] |> copy    # atmospheric water vapor

    # ── Control CO₂ level ───────────────────────────────────────
    CO2_ctrl = cfg.co2_concentration

    if cfg.experiment == :a1b_scenario
        CO2_ctrl = 298.0f0  # A1B scenario baseline
    elseif cfg.experiment in (:a1b_enhanced, :rcp26, :rcp45, :rcp60, :rcp85, :custom_co2,
                               :ssp119, :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2)
        CO2_ctrl = 280.0f0  # IPCC scenarios baseline
    end

    if !cfg.log_co2_dmc
        CO2_ctrl = 0.0f0  # 0 CO2 for deconstruction experiments
    end

    return (Ts_ini=Ts_ini, Ta_ini=Ta_ini, To_ini=To_ini,
        q_ini=q_ini, CO2_ctrl=CO2_ctrl)
end

"""
    qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields, state, timestate, cfg, ws, time_flux; ws_a=ws, ws_q=ws)

Runs `time_flux` years of `tendencies!` to derive the ocean/atmosphere flux
corrections (`fields.TF_correct`/`qF_correct`/`ToF_correct`) that make the
control climate match observed climatology. Mutates `Ts`/`Ta`/`q`/`To` in
place as it integrates. `ws_a`/`ws_q` are forwarded to [`tendencies!`](@ref)
"""
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace, time_flux;
    ws_a::CirculationWorkspace=ws, ws_q::CirculationWorkspace=ws)
    cap_surf = fields.cap_surf
    for it in 1:(time_flux*ndt_days*ndays_yr)
        timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
        timestate.ityr = mod(it - 1, nstep_yr) + 1
        ityr = timestate.ityr

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, fields, state, ws, timestate, cfg; ws_a=ws_a, ws_q=ws_q)

        # Views into climatology & correction fields
        Tc = @view fields.Tclim[:, :, ityr]
        Toc = @view fields.Toclim[:, :, ityr]
        qc = @view fields.qclim[:, :, ityr]
        TFc = @view fields.TF_correct[:, :, ityr]
        ToFc = @view fields.ToF_correct[:, :, ityr]
        qFc = @view fields.qF_correct[:, :, ityr]

        # Surface/air temperature, deep ocean, and humidity update.
        Ts0_buf = ws.Ts0_buf; Ta0_buf = ws.Ta0_buf; To0_buf = ws.To0_buf; q0_buf = ws.q0_buf
        dT_ocean = tend.dT_ocean; SW = tend.SW; LW_surf = tend.LW_surf; LW_down = tend.LW_down
        Q_lat = tend.Q_lat; Q_sens = tend.Q_sens; dTa_crcl = tend.dTa_crcl
        LW_up = tend.LW_up; em = tend.em; Q_lat_air = tend.Q_lat_air
        dTo = tend.dTo; dq_crcl = tend.dq_crcl; dq_eva = tend.dq_eva; dq_rain = tend.dq_rain

        @turbo for j in 1:ydim
            for i in 1:xdim
                ts0 = Ts[i, j] + dT_ocean[i, j] + Δt * (SW[i, j] + LW_surf[i, j] - LW_down[i, j] +
                    Q_lat[i, j] + Q_sens[i, j]) / cap_surf[i, j]
                tfc = (Tc[i, j] - ts0) * cap_surf[i, j] / Δt
                ts0 = ts0 + tfc * Δt / cap_surf[i, j]
                TFc[i, j] = tfc
                Ts0_buf[i, j] = ts0

                Ta0_buf[i, j] = Ta[i, j] + dTa_crcl[i, j] + ΔT_AIR_FACTOR * (
                    LW_up[i, j] + LW_down[i, j] - em[i, j] * LW_surf[i, j] +
                    Q_lat_air[i, j] - Q_sens[i, j])

                to0 = To[i, j] + dTo[i, j]
                tofc = Toc[i, j] - to0
                to0 = to0 + tofc
                ToFc[i, j] = tofc
                To0_buf[i, j] = to0

                q0 = q[i, j] + dq_crcl[i, j] + Δt * (dq_eva[i, j] + dq_rain[i, j])
                qfc = qc[i, j] - q0
                q0 = q0 + qfc
                qFc[i, j] = qfc
                q0_buf[i, j] = q0
            end
        end

        # Sea ice (updates cap_surf in place)
        seaice!(ws.Ts0_buf, fields, timestate, cfg)

        # Diagnostics
        surf = SurfaceState(ws.Ts0_buf, ws.Ta0_buf, ws.To0_buf, ws.q0_buf)
        diagnostics!(it, 0.0, CO2_ctrl, surf, tend, fields, state, timestate)

        # Advance state
        @. Ts = ws.Ts0_buf
        @. Ta = ws.Ta0_buf
        @. q = ws.q0_buf
        @. To = ws.To0_buf
    end
    return nothing
end

const _SOLAR_SWAP_FORCING_TYPE = Dict(
    :paleo_231kyr => :paleo,
    :paleo_solar_modern_co2 => :paleo,
    :obliquity => :obliquity,
    :eccentricity => :eccentricity,
)

const _CO2_SCENARIO_SYMBOLS = (:rcp26, :rcp45, :rcp60, :ssp119, :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2)

const _CO2_SCENARIO_KEY = Dict(:historical_co2 => :hist, :rcp60 => :rcp6)

"""
    greb_model!(run::RunSpec, cfg::PhysicsConfig; jld2_dir="", fields=ClimateFields(),
                allow_uninitialized=false)

Run a GREB flux-correction spin-up (`run.flux` years), control run
(`run.ctrl` years), and scenario run (`run.scnr` years) for `cfg`.

`fields` holds the loaded climatology/grid/flux-correction state (see
[`ClimateFields`](@ref), built by [`load_greb_jld2!`](@ref)). Pass the same
`fields` instance across multiple calls to reuse already-loaded climatology
instead of reloading it - that's the only case where `co2_part`/`sw_solar`
mutations from one run could otherwise leak into the next; this function
resets/restores them per-run regardless.
"""
function greb_model!(run::RunSpec, cfg::PhysicsConfig;
    jld2_dir::AbstractString="", fields::ClimateFields=ClimateFields(),
    allow_uninitialized::Bool=false)
    if !fields.loaded && !allow_uninitialized
        error("""
              greb_model! was given an uninitialized ClimateFields (all-zero climatology).
              
              Load the data and pass it through:
                  fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
                  greb_model!(run, cfg; jld2_dir=jld2_dir, fields=fields)

              If a data-free run is intended (precompilation, config/scenario-plumbing
              tests), opt in explicitly with `allow_uninitialized=true`.
              """)
    end
    time_flux, time_ctrl, time_scnr = run.flux, run.ctrl, run.scnr
    sw_solar_backup = copy(fields.sw_solar)
    try

    state = ModelState()

    # ── 1. Initialisation ───────────────────────────────────────
    if cfg.experiment == :rcp85
        load_cc_anomaly_jld2!(jld2_dir, fields, cfg)
    elseif cfg.experiment == :elnino
        load_enso_anomaly_jld2!(jld2_dir, fields, cfg, :elnino)
    elseif cfg.experiment == :lanina
        load_enso_anomaly_jld2!(jld2_dir, fields, cfg, :lanina)
    end

    ini = init_model!(cfg, fields)
    Ts_ini = ini.Ts_ini;
    Ta_ini = ini.Ta_ini
    To_ini = ini.To_ini;
    q_ini = ini.q_ini
    CO2_ctrl = ini.CO2_ctrl

    # Determine experiment type from cfg
    is_orbital_exp = cfg.experiment in (:obliquity, :eccentricity, :earth_sun_distance)
    is_forced_boundary = cfg.experiment in (:rcp85, :elnino, :lanina)
    is_sst_plus1 = cfg.experiment == :sst_plus1
    is_historical_exp = cfg.experiment == :historical_co2

    # Workspace and accumulator. `ws_a`/`ws_q` are separate `circulation!`
    # scratch spaces so the Ta/q circulation calls inside `tendencies!` can
    # run concurrently on `Threads.nthreads() > 1`
    ws = CirculationWorkspace()
    ws_a = CirculationWorkspace()
    ws_q = CirculationWorkspace()
    acc = MonthlyAccumulator()

    # Initialize time state
    timestate = TimeState(1, 1)

    # ── 2. Flux-correction spin-up ──────────────────────────────
    if !cfg.log_topo_drsp && cfg.log_qflux_dmc
        println("% loading flux correction fields...")
        load_flux_corrections_jld2!(jld2_dir, fields)
    elseif cfg.log_topo_drsp || cfg.log_qflux_dmc
        println("% flux correction  CO2 = ", CO2_ctrl)
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, fields, state, timestate, cfg, ws, time_flux;
            ws_a=ws_a, ws_q=ws_q)
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
    state.sw_solar_forcing = 1.0f0
    mon = 1;
    year = 1970;
    irec = 0

    ctrl_output = MonthlyRecord[]
    sizehint!(ctrl_output, time_ctrl * 12)  # Pre-allocate for all months
    timestate = TimeState(1, 1)  # Initialize time state

    for it in 1:(time_ctrl*nstep_yr)
        (mon, irec) = time_loop!(it, year, CO2_ctrl, mon, irec,
            Ts, Ta, q, To, ctrl_output, fields, state, ws, acc, timestate, cfg;
            ws_a=ws_a, ws_q=ws_q)
        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

    # ── Build ice climatology from control output ───────────────
    ice_forcing = compute_annual_ice_climatology(ctrl_output)

    # ── 4. Scenario run ─────────────────────────────────────────
    println("SCENARIO: ", cfg.experiment, "  time = ", time_scnr, " yr")

    # Paleo/orbital experiments: swap in the alternate solar-forcing table
    if haskey(_SOLAR_SWAP_FORCING_TYPE, cfg.experiment)
        forcing_type = _SOLAR_SWAP_FORCING_TYPE[cfg.experiment]
        println("% loading alternate solar forcing ($forcing_type)...")
        fields.sw_solar .= load_solar_forcing_jld2(jld2_dir, forcing_type, cfg.orbital_index)
    end

    # IPCC RCP/SSP/historical experiments: load the year→CO2 lookup table
    if cfg.experiment in _CO2_SCENARIO_SYMBOLS
        println("% loading IPCC CO2 scenario ($(cfg.experiment))...")
        scenario_key = get(_CO2_SCENARIO_KEY, cfg.experiment, cfg.experiment)
        cfg.co2_scenario = load_co2_scenario_jld2(jld2_dir, scenario_key)
    elseif cfg.experiment == :custom_co2
        isempty(cfg.custom_co2_path) &&
            error("PhysicsConfig.custom_co2_path must be set for the :custom_co2 experiment")
        println("% loading custom CO2 scenario ($(cfg.custom_co2_path))...")
        cfg.co2_scenario = load_custom_co2_scenario(cfg.custom_co2_path)
    end

    # Reset state to initial conditions
    Ts .= Ts_ini;
    Ta .= Ta_ini
    q .= q_ini;
    To .= To_ini
    year = is_orbital_exp ? 1 : (is_historical_exp ? 1850 : 1950)
    CO2 = 340.0f0;
    mon = 1;
    irec = 0

    state.sw_solar_forcing = 1.0f0
    reset!(acc)  # Use accumulator reset

    scnr_output = MonthlyRecord[]
    if time_scnr > 0
        sizehint!(scnr_output, time_scnr * 12)
    end

    for it in 1:(time_scnr*nstep_yr)
        # Obtain forcing (CO2 and solar multiplier)
        forcing_result = forcing(it, year, cfg, fields, ice_forcing; nstep_yr=nstep_yr)
        CO2 = forcing_result.CO2
        state.sw_solar_forcing = forcing_result.sw_solar_forcing

        # Forced‑boundary experiments: overwrite Ts with climatology
        if is_forced_boundary
            ityr_now = mod(it - 1, nstep_yr) + 1
            Ts .= @view fields.Tclim[:, :, ityr_now]
        end

        # SST+1 K experiment
        if is_sst_plus1
            CO2 = CO2_ctrl
            ityr_now = mod(it - 1, nstep_yr) + 1
            @. Ts = ifelse(fields.z_topo < 0.0f0, fields.Tclim[:, :, ityr_now] + 1.0f0, Ts)
        end

        (mon, irec) = time_loop!(it, year, CO2, mon, irec,
            Ts, Ta, q, To, scnr_output, fields, state, ws, acc, timestate, cfg;
            ws_a=ws_a, ws_q=ws_q)

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

    finally
        fields.sw_solar .= sw_solar_backup
    end
end

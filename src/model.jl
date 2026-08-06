# ── notebook cell d404043f-8080-4262-9ab6-d9bb13eee504  (orig lines 1200-1304) ──
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
    fields.co2_part .= 1.0

    # Four of the six regional_co2_* masks depend only on fixed latitude-row
    # ranges — static for the whole run — so compute them once here instead
    # of every timestep in forcing() (IMPROVEMENTS.md §4.4). The other two,
    # :regional_co2_ocean/:regional_co2_land_ice, depend on the control run's
    # ice climatology (icmn_ctrl, src/model.jl's `ice_forcing`), which doesn't
    # exist yet at this point in the run — those stay computed in forcing().
    if cfg.experiment == :regional_co2_nh
        fields.co2_part[:, 1:24] .= 0.5
    elseif cfg.experiment == :regional_co2_sh
        fields.co2_part[:, 25:48] .= 0.5
    elseif cfg.experiment == :regional_co2_tropics
        fields.co2_part[:, 1:15] .= 0.5
        fields.co2_part[:, 33:48] .= 0.5
        for i in 4:4:96
            fields.co2_part[i, 33] = 1.0
            fields.co2_part[i, 15] = 1.0
        end
    elseif cfg.experiment == :regional_co2_extratropics
        fields.co2_part[:, 16:32] .= 0.5
        for i in 4:4:96
            fields.co2_part[i, 32] = 1.0
            fields.co2_part[i, 16] = 1.0
        end
    end

    Tclim = fields.Tclim
    z_topo = fields.z_topo

    # ── dTrad: offset between T_atm and radiation temperature ────
    @. fields.dTrad = -0.16 * Tclim - 5.0

    # ── z_ocean: 3× maximum mixed-layer depth over the year ──────
    fields.z_ocean .= 3.0 .* dropdims(maximum(fields.mldclim; dims=3); dims=3)

    # ── Sensitivity experiment overrides ─────────────────────────
    if !cfg.log_topo_drsp
        @. z_topo = min(z_topo, 1.0)       # constant topography
    end

    # Apply cloud deconstruction switch
    if !cfg.log_clouds_dmc
        fields.cldclim .= 0.0  # zero cloud climatology
    end

    # Apply flux correction conditional zeroing (MSCM feature)
    if !cfg.log_topo_drsp && !cfg.log_qflux_dmc
        fields.TF_correct .= 0.0
        fields.qF_correct .= 0.0
        fields.ToF_correct .= 0.0
    end

    # Climatology modifications
    if !cfg.log_hydro_dmc
        fields.qclim .= 0.0  # zero out humidity climatology
    end

    if !cfg.log_clouds_drsp
        fields.cldclim .= 0.7           # constant cloud cover (2xCO2 deconstruction)
    end

    if !cfg.log_humid_drsp
        fields.qclim .= 0.0052          # constant water vapor
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

    # ── Surface heat capacity ────────────────────────────────────
    cap_surf = fields.cap_surf
    mldclim = fields.mldclim
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
    To_ini = fields.Toclim[:, :, nstep_yr] |> copy   # deep ocean temperature
    q_ini = fields.qclim[:, :, nstep_yr] |> copy   # atmospheric water vapor

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
"""
    qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields, state, timestate, cfg, ws, time_flux)

Runs `time_flux` years of `tendencies!` to derive the ocean/atmosphere flux
corrections (`fields.TF_correct`/`qF_correct`/`ToF_correct`) that make the
control climate match observed climatology. Mutates `Ts`/`Ta`/`q`/`To` in
place as it integrates.
"""
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, fields::ClimateFields, state::ModelState, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace, time_flux)
    cap_surf = fields.cap_surf
    for it in 1:(time_flux*ndt_days*ndays_yr)
        timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
        timestate.ityr = mod(it - 1, nstep_yr) + 1
        ityr = timestate.ityr

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, fields, state, ws, timestate, cfg)

        # Views into climatology & correction fields
        Tc = @view fields.Tclim[:, :, ityr]
        Toc = @view fields.Toclim[:, :, ityr]
        qc = @view fields.qclim[:, :, ityr]
        TFc = @view fields.TF_correct[:, :, ityr]
        ToFc = @view fields.ToF_correct[:, :, ityr]
        qFc = @view fields.qF_correct[:, :, ityr]

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
        # No flux-correction term here (unlike Ts/To/q below) — this is
        # intentional, not a gap: Fortran's own qflux_correction never
        # nudges Ta toward a climatology either (greb.model.mscm.f90:566,
        # `Ta0 = Ta1 + dTa + dTa_crcl`), and there is no `TaF_correct`
        # array anywhere in the Fortran source. Verified directly, not
        # assumed — see the regression test asserting this stays true.
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
        seaice!(ws.Ts0_buf, fields, timestate, cfg)

        # Diagnostics
        surf = SurfaceState(ws.Ts0_buf, ws.Ta0_buf, ws.To0_buf, ws.q0_buf)
        diagnostics!(it, 0.0, CO2_ctrl, surf, tend, fields, state, timestate)

        # Advance state (broadcast – same as `.=`)
        @. Ts = ws.Ts0_buf
        @. Ta = ws.Ta0_buf
        @. q = ws.q0_buf
        @. To = ws.To0_buf
    end
    return nothing
end

# ── notebook cell 92c3bd68-bd07-4381-9c04-e6611650cd1e  (orig lines 2929-3043) ──
# Experiments that swap in a full alternate (ydim, nstep_yr) solar-forcing
# table, mirroring Fortran's `sw_solar = sw_solar_scnr` (log_exp 30/31/35/36).
const _SOLAR_SWAP_FORCING_TYPE = Dict(
    :paleo_231kyr => :paleo,
    :paleo_solar_modern_co2 => :paleo,
    :obliquity => :obliquity,
    :eccentricity => :eccentricity,
)

"""
    greb_model!(run::RunSpec, cfg::PhysicsConfig; jld2_dir="", fields=ClimateFields())

Run a GREB flux-correction spin-up (`run.flux` years), control run
(`run.ctrl` years), and scenario run (`run.scnr` years) for `cfg`.

`fields` holds the loaded climatology/grid/flux-correction state (see
[`ClimateFields`](@ref), built by [`load_greb_jld2!`](@ref)); it defaults to
a fresh all-zero instance so callers that never load real data (tests,
quick structural runs) don't need to pass anything. Pass the same `fields`
instance across multiple calls to reuse already-loaded climatology instead
of reloading it — that's the only case where `co2_part`/`sw_solar`
mutations from one run could otherwise leak into the next; this function
resets/restores them per-run regardless.
"""
function greb_model!(run::RunSpec, cfg::PhysicsConfig;
    jld2_dir::AbstractString="", fields::ClimateFields=ClimateFields())
    time_flux, time_ctrl, time_scnr = run.flux, run.ctrl, run.scnr
    # Back up sw_solar so a paleo/orbital run's swapped table never leaks
    # into a later run against the same (possibly reused) `fields` instance.
    sw_solar_backup = copy(fields.sw_solar)
    try

    state = ModelState()

    # ── 1. Initialisation ───────────────────────────────────────
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

    # Workspace and accumulator
    ws = CirculationWorkspace()
    acc = MonthlyAccumulator()

    # Initialize time state
    timestate = TimeState(1, 1)

    # ── 2. Flux-correction spin-up ──────────────────────────────
    if cfg.log_topo_drsp || cfg.log_qflux_dmc
        if !cfg.log_topo_drsp && cfg.log_qflux_dmc
            println("% loading flux correction fields...")
            load_flux_corrections_jld2!(jld2_dir, fields)
        end
        println("% flux correction  CO2 = ", CO2_ctrl)
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, fields, state, timestate, cfg, ws, time_flux)
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
    state.sw_solar_forcing = 1.0
    mon = 1;
    year = 1970;
    irec = 0

    ctrl_output = MonthlyRecord[]
    sizehint!(ctrl_output, time_ctrl * 12)  # Pre-allocate for all months
    timestate = TimeState(1, 1)  # Initialize time state

    for it in 1:(time_ctrl*nstep_yr)
        (mon, irec) = time_loop!(it, year, CO2_ctrl, mon, irec,
            Ts, Ta, q, To, ctrl_output, fields, state, ws, acc, timestate, cfg)
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

    # Paleo/orbital experiments: swap in the alternate solar-forcing table
    # (Fortran: `sw_solar = sw_solar_scnr` at the start of the scenario run).
    # :modern_solar_paleo_co2 intentionally does NOT swap solar — matches
    # Fortran's log_exp==32, which only changes CO2.
    if haskey(_SOLAR_SWAP_FORCING_TYPE, cfg.experiment)
        forcing_type = _SOLAR_SWAP_FORCING_TYPE[cfg.experiment]
        println("% loading alternate solar forcing ($forcing_type)...")
        fields.sw_solar .= load_solar_forcing_jld2(jld2_dir, forcing_type, cfg.orbital_index)
    end

    # Reset state to initial conditions
    Ts .= Ts_ini;
    Ta .= Ta_ini
    q .= q_ini;
    To .= To_ini
    year = is_orbital_exp ? 1 : 1950
    CO2 = 340.0;
    mon = 1;
    irec = 0

    state.sw_solar_forcing = 1.0
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
            @. Ts = ifelse(fields.z_topo < 0.0, fields.Tclim[:, :, ityr_now] + 1.0, Ts)
        end

        (mon, irec) = time_loop!(it, year, CO2, mon, irec,
            Ts, Ta, q, To, scnr_output, fields, state, ws, acc, timestate, cfg)

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

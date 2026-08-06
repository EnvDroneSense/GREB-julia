"""
    diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, fields, state, timestate)

Accumulates the current timestep into `state`'s annual-mean buffers; at the
last timestep of the year, averages them, prints the annual summary line
(global mean + two sample points), and resets the accumulators for the next
year.
"""
# ── notebook cell cc0e682c-8767-498e-8178-2de4e796b3a8  (orig lines 2355-2392) ──
function diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, fields::ClimateFields, state::ModelState, timestate)
    # Accumulate
    state.Tsmn .+= Ts0;
    state.Tamn .+= Ta0;
    state.Tomn .+= To0
    state.qmn .+= q0;
    state.amn .+= albedo
    state.swmn .+= sw;
    state.lwmn .+= lw_surf
    state.qlatmn .+= q_lat;
    state.qsensmn .+= q_sens
    state.ftmn .+= @view fields.TF_correct[:, :, timestate.ityr]
    state.fqmn .+= @view fields.qF_correct[:, :, timestate.ityr]

    if timestate.ityr == nstep_yr
        # Compute annual means
        n = nstep_yr
        state.Tsmn ./= n;
        state.Tamn ./= n;
        state.Tomn ./= n
        state.qmn ./= n;
        state.amn ./= n
        state.swmn ./= n;
        state.lwmn ./= n
        state.qlatmn ./= n;
        state.qsensmn ./= n
        state.ftmn ./= n;
        state.fqmn ./= n

        # Global mean and sample points (°C)
        global_mean = sum(state.Tsmn) / (xdim * ydim) - 273.15
        point1 = state.Tsmn[48, 27] - 273.15   # Tropical Pacific
        point2 = state.Tsmn[16, 38] - 273.15   # Hamburg/North Europe

        println(year, "  ", round(global_mean, digits=2),
            "  ", round(point1, digits=2),
            "  ", round(point2, digits=2))

        # Reset accumulators
        fill!(state.Tsmn, 0.0);
        fill!(state.Tamn, 0.0);
        fill!(state.Tomn, 0.0)
        fill!(state.qmn, 0.0);
        fill!(state.amn, 0.0)
        fill!(state.swmn, 0.0);
        fill!(state.lwmn, 0.0)
        fill!(state.qlatmn, 0.0);
        fill!(state.qsensmn, 0.0)
        fill!(state.ftmn, 0.0);
        fill!(state.fqmn, 0.0)
    end
    return nothing
end

"""
    output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens, output_buf, acc, timestate)

Accumulates the current timestep into `acc`; on the last timestep of `mon`,
pushes a monthly-mean [`MonthlyRecord`](@ref) onto `output_buf`, resets `acc`,
and advances to the next month. Returns `(irec, mon)`.
"""
# ── notebook cell dfdde9f1-b226-4af0-9ac2-36f1b01622fa  (orig lines 2405-2437) ──
function output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens,
    output_buf::Vector{MonthlyRecord}, acc::MonthlyAccumulator, timestate)
    # ----- SAFETY: clamp month to 1..12 -----
    mon = clamp(mon, 1, 12)

    accumulate!(acc, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)

    # ----- Check end of month -----
    if timestate.jday == jday_mon_cumsum[mon] && (it % ndt_days == 0)
        ndm = cjday_mon[mon] * ndt_days
        irec += 1
        push!(output_buf, (
            Ts=copy(acc.Tmm ./ ndm),
            Ta=copy(acc.Tamm ./ ndm),
            To=copy(acc.Tomm ./ ndm),
            q=copy(acc.qmm ./ ndm),
            albedo=copy(acc.apmm ./ ndm),
            ice=copy(acc.icemm ./ ndm),
            precip=copy(acc.precipmm ./ ndm),
            evap=copy(acc.evapmm ./ ndm),
            qcrcl=copy(acc.qcrclmm ./ ndm),
            sw=copy(acc.swmm ./ ndm),
            lw=copy(acc.lwmm ./ ndm),
            qlat=copy(acc.qlatmm ./ ndm),
            qsens=copy(acc.qsensmm ./ ndm)
        ))
        reset!(acc)
        mon = mon == 12 ? 1 : mon + 1
    end
    return (irec=irec, mon=mon)
end

"""
    time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf, fields, state, ws, acc, timestate, cfg)

One full model timestep: computes [`tendencies!`](@ref), integrates
`Ts`/`Ta`/`To`/`q` forward with flux corrections applied, runs
[`seaice!`](@ref), then dispatches to [`output!`](@ref) and
[`diagnostics!`](@ref). Returns `(mon, irec)`.
"""
# ── notebook cell 33fa7b1f-938b-481c-bf25-eca8d7fb33a7  (orig lines 2438-2500) ──
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
    fields::ClimateFields, state::ModelState, ws::CirculationWorkspace, acc::MonthlyAccumulator,
    timestate, cfg::PhysicsConfig)
    # Calendar lookup
    cal = it <= max_timesteps ? calendar_lookup[it] : (
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    )
    timestate.jday = cal.day
    timestate.ityr = cal.step
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)

    # Correction views
    TF_corr = @view fields.TF_correct[:, :, ityr]
    qF_corr = @view fields.qF_correct[:, :, ityr]
    ToF_corr = @view fields.ToF_correct[:, :, ityr]
    cap_surf = fields.cap_surf
    wz_vapor = fields.wz_vapor

    # Surface temperature
    @. Ts = Ts + tend.dT_ocean + Δt * (tend.SW + tend.LW_surf - tend.LW_down +
                                       tend.Q_lat + tend.Q_sens + TF_corr) / cap_surf

    # Air temperature
    @. Ta = Ta + tend.dTa_crcl + Δt * (tend.LW_up + tend.LW_down - tend.em *
                                                                   tend.LW_surf + tend.Q_lat_air - tend.Q_sens) / cap_air

    # Clamps
    @. Ts = max(Ts, min_T_K)
    @. Ta = max(Ta, min_T_K)

    # Deep ocean
    @. To = To + tend.dTo + ToF_corr

    # Humidity: Fortran does two threshold-triggered replacements, not a
    # standard clamp (greb.model.mscm.f90:481-483) — `dq<=-q1` replaces with
    # -0.9*q1 (leaving dq UNCHANGED for the range between -q1 and -0.9*q1,
    # unlike `clamp` which would pull anything below -0.9*q1 up to it), and
    # `dq>0.020` replaces with 0.020. Also, when log_hydro_dmc==0 Fortran
    # zeroes the ENTIRE dq afterward (greb.model.mscm.f90:486) — not just the
    # eva/rain terms (already zero from hydro!'s own early return) but
    # dq_crcl and qF_correct too, which the old `clamp`-based version never
    # zeroed. tend.dq_eva/dq_rain are already the right values regardless of
    # log_hydro_dmc (hydro! itself returns zero buffers when it's off).
    dq_eva_use = tend.dq_eva
    dq_rain_use = tend.dq_rain
    dq_crcl_use = cfg.log_crcl_dmc ? tend.dq_crcl : ws.crcl
    @. ws.temp_buf = Δt * (dq_eva_use + dq_rain_use) + dq_crcl_use + qF_corr
    @. ws.temp_buf = ifelse(ws.temp_buf <= -q, -min_humidity_change * q, ws.temp_buf)
    @. ws.temp_buf = ifelse(ws.temp_buf > max_humidity_change, max_humidity_change, ws.temp_buf)
    if !cfg.log_hydro_dmc
        fill!(ws.temp_buf, 0.0)
    end
    @. q = q + ws.temp_buf

    # Sea ice heat capacity
    seaice!(Ts, fields, timestate, cfg)

    # Conversion to mm/day (analysis units)
    @. ws.precip_out = (-dq_rain_use) * wz_vapor * conv_factor
    @. ws.evap_out = dq_eva_use * wz_vapor * conv_factor
    @. ws.qcrcl_out = dq_crcl_use


    # Output and diagnostics
    (mon, irec) = output!(it, irec, mon, Ts, Ta, To, q, tend.albedo,
        tend.ice_cover, ws.precip_out, ws.evap_out, ws.qcrcl_out,
        tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens,
        output_buf, acc, timestate)
    diagnostics!(it, year, CO2, Ts, Ta, To, q, tend.albedo,
        tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, fields, state, timestate)

    return (mon=mon, irec=irec)
end

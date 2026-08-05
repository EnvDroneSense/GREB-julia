# ── notebook cell cc0e682c-8767-498e-8178-2de4e796b3a8  (orig lines 2355-2392) ──
function diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate)
    # Accumulate
    Tsmn .+= Ts0;
    Tamn .+= Ta0;
    Tomn .+= To0
    qmn .+= q0;
    amn .+= albedo
    swmn .+= sw;
    lwmn .+= lw_surf
    qlatmn .+= q_lat;
    qsensmn .+= q_sens
    ftmn .+= @view TF_correct[:, :, timestate.ityr]
    fqmn .+= @view qF_correct[:, :, timestate.ityr]

    if timestate.ityr == nstep_yr
        # Compute annual means
        n = nstep_yr
        Tsmn ./= n;
        Tamn ./= n;
        Tomn ./= n
        qmn ./= n;
        amn ./= n
        swmn ./= n;
        lwmn ./= n
        qlatmn ./= n;
        qsensmn ./= n
        ftmn ./= n;
        fqmn ./= n

        # Global mean and sample points (°C)
        global_mean = sum(Tsmn) / (xdim * ydim) - 273.15
        point1 = Tsmn[48, 27] - 273.15   # Tropical Pacific
        point2 = Tsmn[16, 38] - 273.15   # Hamburg/North Europe

        println(year, "  ", round(global_mean, digits=2),
            "  ", round(point1, digits=2),
            "  ", round(point2, digits=2))

        # Reset accumulators
        fill!(Tsmn, 0.0);
        fill!(Tamn, 0.0);
        fill!(Tomn, 0.0)
        fill!(qmn, 0.0);
        fill!(amn, 0.0)
        fill!(swmn, 0.0);
        fill!(lwmn, 0.0)
        fill!(qlatmn, 0.0);
        fill!(qsensmn, 0.0)
        fill!(ftmn, 0.0);
        fill!(fqmn, 0.0)
    end
    return nothing
end

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

# ── notebook cell 33fa7b1f-938b-481c-bf25-eca8d7fb33a7  (orig lines 2438-2500) ──
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
    ws::CirculationWorkspace, acc::MonthlyAccumulator,
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
    tend = tendencies!(CO2, Ts, Ta, To, q, ws, timestate, cfg)

    # Correction views
    TF_corr = @view TF_correct[:, :, ityr]
    qF_corr = @view qF_correct[:, :, ityr]
    ToF_corr = @view ToF_correct[:, :, ityr]

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

    # Humidity (with clamp)
    dq_eva_use = cfg.log_hydro_dmc ? tend.dq_eva : ws.eva
    dq_rain_use = cfg.log_hydro_dmc ? tend.dq_rain : ws.rain
    dq_crcl_use = cfg.log_crcl_dmc ? tend.dq_crcl : ws.crcl
    @. q = q + clamp(Δt * (dq_eva_use + dq_rain_use) + dq_crcl_use + qF_corr,
        -min_humidity_change * q, max_humidity_change)

    # Sea ice heat capacity
    seaice!(Ts, timestate, cfg)

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
        tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

    return (mon=mon, irec=irec)
end

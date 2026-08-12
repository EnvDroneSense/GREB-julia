"""
    diagnostics!(it, year, CO2, surf::SurfaceState, tend, fields, state, timestate)

Accumulates the current timestep into `state`'s annual-mean buffers; at the
last timestep of the year, averages them, prints the annual summary line
(global mean + two sample points), and resets the accumulators for the next
year. `tend` is the `NamedTuple` [`tendencies!`](@ref) returns.
"""
function diagnostics!(it, year, CO2, surf::SurfaceState, tend, fields::ClimateFields, state::ModelState, timestate)
    # Accumulate
    Tsmn = state.Tsmn; Tamn = state.Tamn; Tomn = state.Tomn; qmn = state.qmn
    amn = state.amn; swmn = state.swmn; lwmn = state.lwmn
    qlatmn = state.qlatmn; qsensmn = state.qsensmn
    ftmn = state.ftmn; fqmn = state.fqmn
    Ts = surf.Ts; Ta = surf.Ta; To = surf.To; q = surf.q
    albedo = tend.albedo; SW = tend.SW; LW_surf = tend.LW_surf
    Q_lat = tend.Q_lat; Q_sens = tend.Q_sens
    TF_correct = fields.TF_correct; qF_correct = fields.qF_correct
    ityr = timestate.ityr

    @turbo for j in 1:ydim
        for i in 1:xdim
            Tsmn[i, j] += Ts[i, j]
            Tamn[i, j] += Ta[i, j]
            Tomn[i, j] += To[i, j]
            qmn[i, j] += q[i, j]
            amn[i, j] += albedo[i, j]
            swmn[i, j] += SW[i, j]
            lwmn[i, j] += LW_surf[i, j]
            qlatmn[i, j] += Q_lat[i, j]
            qsensmn[i, j] += Q_sens[i, j]
            ftmn[i, j] += TF_correct[i, j, ityr]
            fqmn[i, j] += qF_correct[i, j, ityr]
        end
    end

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
        global_mean = sum(state.Tsmn[i, j] * dxlat_grid[j] for i in 1:xdim, j in 1:ydim) /
                      (xdim * sum(dxlat_grid)) - 273.15
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
    output!(it, irec, mon, surf::SurfaceState, tend, ws, output_buf, acc, timestate)

Accumulates the current timestep into `acc`; on the last timestep of `mon`,
pushes a monthly-mean [`MonthlyRecord`](@ref) onto `output_buf`, resets `acc`,
and advances to the next month. Returns `(mon, irec)`. `tend` is the
`NamedTuple` [`tendencies!`](@ref) returns; `ws.precip_out`/`evap_out`/
`qcrcl_out` hold this step's converted precipitation/evaporation/moisture-
circulation output.
"""
function output!(it, irec, mon, surf::SurfaceState, tend, ws::CirculationWorkspace,
    output_buf::Vector{MonthlyRecord}, acc::MonthlyAccumulator, timestate)
    mon = clamp(mon, 1, 12)

    accumulate!(acc, surf.Ts, surf.Ta, surf.To, surf.q, tend.albedo, tend.ice_cover,
        ws.precip_out, ws.evap_out, ws.qcrcl_out, tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens)

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
    return (mon=mon, irec=irec)
end

"""
    time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf, fields, state, ws, acc, timestate, cfg; ws_a=ws, ws_q=ws)

One full model timestep: computes [`tendencies!`](@ref), integrates
`Ts`/`Ta`/`To`/`q` forward with flux corrections applied, runs
[`seaice!`](@ref), then dispatches to [`output!`](@ref) and
[`diagnostics!`](@ref). Returns `(mon, irec)`. `ws_a`/`ws_q` are forwarded to
[`tendencies!`](@ref) - see its docstring for the opt-in threading they
enable.
"""
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
    fields::ClimateFields, state::ModelState, ws::CirculationWorkspace, acc::MonthlyAccumulator,
    timestate, cfg::PhysicsConfig; ws_a::CirculationWorkspace=ws, ws_q::CirculationWorkspace=ws)
    # Calendar lookup
    cal = it <= max_timesteps ? calendar_lookup[it] : (
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    )
    timestate.jday = cal.day
    timestate.ityr = cal.step
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg; ws_a=ws_a, ws_q=ws_q)

    # Correction views
    TF_corr = @view fields.TF_correct[:, :, ityr]
    qF_corr = @view fields.qF_correct[:, :, ityr]
    ToF_corr = @view fields.ToF_correct[:, :, ityr]
    cap_surf = fields.cap_surf
    wz_vapor = fields.wz_vapor

    # Humidity tendency buffer selection
    dq_eva_use = tend.dq_eva
    dq_rain_use = tend.dq_rain
    dq_crcl_use = cfg.log_crcl_dmc ? tend.dq_crcl : ws.crcl
    hydro_on = cfg.log_hydro_dmc ? 1.0 : 0.0

    SW = tend.SW; LW_surf = tend.LW_surf; LW_down = tend.LW_down
    Q_lat = tend.Q_lat; Q_sens = tend.Q_sens; dTa_crcl = tend.dTa_crcl
    LW_up = tend.LW_up; em = tend.em; Q_lat_air = tend.Q_lat_air
    dTo = tend.dTo; dT_ocean = tend.dT_ocean
    temp_buf = ws.temp_buf; precip_out = ws.precip_out
    evap_out = ws.evap_out; qcrcl_out = ws.qcrcl_out

    # Surface/air temperature, deep ocean, and humidity update
    @turbo for j in 1:ydim
        for i in 1:xdim
            Ts[i, j] = Ts[i, j] + dT_ocean[i, j] + Δt * (SW[i, j] + LW_surf[i, j] - LW_down[i, j] +
                Q_lat[i, j] + Q_sens[i, j] + TF_corr[i, j]) / cap_surf[i, j]
            Ta[i, j] = Ta[i, j] + dTa_crcl[i, j] + Δt * (LW_up[i, j] + LW_down[i, j] -
                em[i, j] * LW_surf[i, j] + Q_lat_air[i, j] - Q_sens[i, j]) / cap_air

            Ts[i, j] = max(Ts[i, j], min_T_K)
            Ta[i, j] = max(Ta[i, j], min_T_K)

            To[i, j] = To[i, j] + dTo[i, j] + ToF_corr[i, j]

            tb = Δt * (dq_eva_use[i, j] + dq_rain_use[i, j]) + dq_crcl_use[i, j] + qF_corr[i, j]
            tb = ifelse(tb <= -q[i, j], -min_humidity_change * q[i, j], tb)
            tb = ifelse(tb > max_humidity_change, max_humidity_change, tb)
            tb = hydro_on * tb
            temp_buf[i, j] = tb
            q[i, j] = q[i, j] + tb

            precip_out[i, j] = (-dq_rain_use[i, j]) * wz_vapor[i, j] * conv_factor
            evap_out[i, j] = dq_eva_use[i, j] * wz_vapor[i, j] * conv_factor
            qcrcl_out[i, j] = dq_crcl_use[i, j]
        end
    end

    # Sea ice heat capacity
    seaice!(Ts, fields, timestate, cfg)


    # Output and diagnostics
    surf = SurfaceState(Ts, Ta, To, q)
    (mon, irec) = output!(it, irec, mon, surf, tend, ws, output_buf, acc, timestate)
    diagnostics!(it, year, CO2, surf, tend, fields, state, timestate)

    return (mon=mon, irec=irec)
end

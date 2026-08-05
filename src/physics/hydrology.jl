# ── notebook cell 606032a2-b2ca-4fd8-9930-afd83aecec7a  (orig lines 1495-1608) ──
function hydro!(Ts, q, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    c_q = cfg.c_q
    c_rq = cfg.c_rq
    c_omega = cfg.c_omega
    c_omegastd = cfg.c_omegastd

    fill!(ws.Q_lat_buf, 0.0)
    fill!(ws.Q_lat_air_buf, 0.0)
    fill!(ws.dq_eva_buf, 0.0)
    fill!(ws.dq_rain_buf, 0.0)

    if !cfg.log_atmos_dmc || !cfg.log_hydro_dmc || !cfg.log_hydro_drsp
        return (Q_lat=ws.Q_lat_buf, Q_lat_air=ws.Q_lat_air_buf,
            dq_eva=ws.dq_eva_buf, dq_rain=ws.dq_rain_buf)
    end

    u = @view uclim[:, :, timestate.ityr]
    v = @view vclim[:, :, timestate.ityr]
    swet = @view swetclim[:, :, timestate.ityr]
    omega = @view omegaclim[:, :, timestate.ityr]
    omegastd = @view omegastdclim[:, :, timestate.ityr]

    const_factor1 = 3.75e-3
    const_factor2 = 17.08085
    const_factor3 = 234.175
    gust_land = 4.0
    gust_ocean = 9.0
    cE_land = 0.25 * ce
    cE_ocean = 0.58 * ce
    const_latent = cq_latent * ρ_air * ce

    # Saturation humidity
    @turbo for j in 1:ydim
        for i in 1:xdim
            T = Ts[i, j] - 273.15
            ws.qs[i, j] = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j]
            ws.qs[i, j] = max(ws.qs[i, j], 1e-8)
        end
    end

    # Evaporation
    if cfg.log_eva == -1
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i, j];
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i, j] = (q[i, j] - ws.qs[i, j]) * wind * const_latent * swet[i, j]
            end
        end
    elseif cfg.log_eva == 0
        ws_view = @view wsclim[:, :, timestate.ityr]
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.Tskin[i, j] = ifelse(z_topo[i, j] > 0.0, Ts[i, j] + 5.0, Ts[i, j] + 1.0)
                ws.Tskin[i, j] = ifelse(ws.Tskin[i, j] < 200.0, 200.0, ws.Tskin[i, j])
                T = ws.Tskin[i, j] - 273.15
                qs_val = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j]
                ws.qs[i, j] = qs_val

                ws.ws_base[i, j] = ws_view[i, j]
                gust = ifelse(z_topo[i, j] > 0.0, 132.25, 29.16)
                wind = sqrt(ws.ws_base[i, j]*ws.ws_base[i, j] + gust)

                ws.cE_buf[i, j] = ifelse(z_topo[i, j] > 0.0, cE_land, cE_ocean)
                ws.Q_lat_buf[i, j] = ws.cE_buf[i, j] * wind * ρ_air * cq_latent * (q[i, j] - qs_val) * swet[i, j]
            end
        end
    else
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i, j];
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i, j] = (q[i, j] - ws.qs[i, j]) * wind * const_latent * swet[i, j]
            end
        end
    end

    # Precipitation - use ws.rq buffer
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.rq[i, j] = q[i, j] / max(ws.qs[i, j], 1e-8)
            ws.dq_rain_buf[i, j] = (c_q + c_rq * ws.rq[i, j] + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
        end
    end

    # Apply rain limit
    if cfg.log_rain == 1
        limit_val = -0.0015 / (wz_vapor[1, 1] * r_qviwv * 86400.0)
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.dq_rain_buf[i, j] = ifelse(ws.dq_rain_buf[i, j] >= limit_val, limit_val, ws.dq_rain_buf[i, j])
            end
        end
    end

    # Water vapor tendencies
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.dq_eva_buf[i, j] = -ws.Q_lat_buf[i, j] / cq_latent / r_qviwv
            min_dq = -0.9 * q[i, j] / Δt
            ws.dq_rain_buf[i, j] = ifelse(ws.dq_rain_buf[i, j] < min_dq, min_dq, ws.dq_rain_buf[i, j])
            ws.Q_lat_air_buf[i, j] = -ws.dq_rain_buf[i, j] * cq_latent * r_qviwv
        end
    end

    return (Q_lat=ws.Q_lat_buf,
        Q_lat_air=ws.Q_lat_air_buf,
        dq_eva=ws.dq_eva_buf,
        dq_rain=ws.dq_rain_buf)
end

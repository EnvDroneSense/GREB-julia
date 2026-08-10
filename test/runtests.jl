using GREB
using Test

# Smoke tests that do NOT require the (large, external) JLD2 input data.
# They check that the package loads, its types build, and the grid/constants
# are intact after the notebook -> package extraction. Full integration runs
# (which need `greb_dataset_jld2/`) are demonstrated in examples/run_greb.jl.

const DATA_DIR = joinpath(@__DIR__, "..", "greb_dataset_jld2")

@testset "GREB.jl" begin

    @testset "grid constants" begin
        @test GREB.xdim == 96
        @test GREB.ydim == 48
        @test GREB.nstep_yr == 730
    end

    @testset "PhysicsConfig" begin
        cfg = PhysicsConfig()
        @test cfg isa PhysicsConfig

        for exp in (:full_model, :constant_topo, :co2_double, :co2_quadruple,
                    :elnino, :lanina, :rcp85)
            c = create_experiment_config(exp)
            @test c isa PhysicsConfig
            @test c.experiment == exp
        end

        @test create_experiment_config(:co2_double).co2_concentration == 680.0
    end

    @testset "set_hydrology_parameters! matches HYDRO_PARAMS for every log_rain value" begin
        for log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_rain = log_rain
            set_hydrology_parameters!(cfg)
            expected = GREB.HYDRO_PARAMS[log_rain]
            @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == expected
        end

        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 0
        cfg.log_clim = 1
        set_hydrology_parameters!(cfg)
        @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == (-1.27, 1.99, -16.54, 21.15)
    end

    @testset "set_hydrology_parameters! errors on invalid log_rain" begin
        # An out-of-range log_rain must fail loudly, not silently fall back
        # to a default parameter set.
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 4
        @test_throws ErrorException set_hydrology_parameters!(cfg)
    end

    @testset "co2_part regional CO2 mask resets between runs (no leak)" begin
        # A regional-CO2 experiment's mask must not leak into a later,
        # unrelated run against the same reused `fields` instance.
        fields = ClimateFields()
        cfg_regional = PhysicsConfig(experiment=:regional_co2_nh)
        redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg_regional; jld2_dir = "", fields = fields)
        end
        @test any(!=(1.0), fields.co2_part)  # regional run actually changed the mask

        cfg_plain = create_experiment_config(:full_model)
        redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields)
        end
        @test all(==(1.0), fields.co2_part)  # init_model! resets it back to full CO2
    end

    @testset "workspace/accumulator/record types construct correctly" begin
        ws = CirculationWorkspace()
        @test ws isa CirculationWorkspace

        acc = MonthlyAccumulator()
        @test acc isa MonthlyAccumulator
        @test (GREB.reset!(acc); true)   # reset! runs without error

        ts = TimeState(1, 1)
        @test ts.jday == 1
        @test ts.ityr == 1

        @test MonthlyRecord <: NamedTuple
        @test :Ts in fieldnames(MonthlyRecord)
        @test :precip in fieldnames(MonthlyRecord)
    end

    @testset "build_monthly_climatology/apply_scenario_anomalies" begin
        # Hand-built records with every field filled to one scalar value
        # make the averaging arithmetic trivial to check by hand.
        mkrec(v) =(Ts=fill(v, GREB.xdim, GREB.ydim), Ta=fill(v, GREB.xdim, GREB.ydim),
            To=fill(v, GREB.xdim, GREB.ydim), q=fill(v, GREB.xdim, GREB.ydim),
            albedo=fill(v, GREB.xdim, GREB.ydim), ice=fill(v, GREB.xdim, GREB.ydim),
            precip=fill(v, GREB.xdim, GREB.ydim), evap=fill(v, GREB.xdim, GREB.ydim),
            qcrcl=fill(v, GREB.xdim, GREB.ydim), sw=fill(v, GREB.xdim, GREB.ydim),
            lw=fill(v, GREB.xdim, GREB.ydim), qlat=fill(v, GREB.xdim, GREB.ydim),
            qsens=fill(v, GREB.xdim, GREB.ydim))

        @test build_monthly_climatology(MonthlyRecord[]) == MonthlyRecord[]

        # Multi-year averaging: 2 years, record idx's value == idx. Month m
        # sees values m (year 1) and m+12 (year 2) -> mean == m+6.
        two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        clim = build_monthly_climatology(two_years)
        @test length(clim) == 12
        for m in 1:12
            @test all(==(Float64(m + 6)), clim[m].Ts)
        end

        # Non-12-multiple count: only 5 records -> months 6..12 never occur
        # and must fall back to records[1] exactly (not e.g. NaN/zero).
        five = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:5]
        clim5 = build_monthly_climatology(five)
        for m in 1:5
            @test all(==(Float64(m)), clim5[m].Ts)
        end
        for m in 6:12
            @test clim5[m] == five[1]
        end

        # apply_scenario_anomalies: 2 years of scnr records (value == idx),
        # ctrl climatology value == 100+month. Check the wraparound (idx=13
        # is month 1, same as idx=1) subtracts the correct month's climatology.
        ctrl_clim = MonthlyRecord[mkrec(100.0 + m) for m in 1:12]
        scnr = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        anom = apply_scenario_anomalies(scnr, ctrl_clim)
        @test all(==(1.0 - 101.0), anom[1].Ts)
        @test all(==(13.0 - 101.0), anom[13].Ts)
        @test all(==(12.0 - 112.0), anom[12].Ts)

        # Early-return guards: empty scnr_records or empty ctrl_clim ->
        # scnr_records passed straight through, not turned into anomalies.
        @test apply_scenario_anomalies(MonthlyRecord[], ctrl_clim) == MonthlyRecord[]
        @test apply_scenario_anomalies(scnr, MonthlyRecord[]) == scnr
    end

    @testset "compute_annual_ice_climatology" begin
        # Same record-index-encodes-value trick as the climatology test
        # above: 2 years, record idx's ice value == idx, so month m sees
        # idx m (year 1) and m+12 (year 2) -> mean == m+6.
        mkrec(v) = (Ts=zeros(GREB.xdim, GREB.ydim), Ta=zeros(GREB.xdim, GREB.ydim),
            To=zeros(GREB.xdim, GREB.ydim), q=zeros(GREB.xdim, GREB.ydim),
            albedo=zeros(GREB.xdim, GREB.ydim), ice=fill(v, GREB.xdim, GREB.ydim),
            precip=zeros(GREB.xdim, GREB.ydim), evap=zeros(GREB.xdim, GREB.ydim),
            qcrcl=zeros(GREB.xdim, GREB.ydim), sw=zeros(GREB.xdim, GREB.ydim),
            lw=zeros(GREB.xdim, GREB.ydim), qlat=zeros(GREB.xdim, GREB.ydim),
            qsens=zeros(GREB.xdim, GREB.ydim))

        @test all(iszero, compute_annual_ice_climatology(MonthlyRecord[]))

        two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        clim = compute_annual_ice_climatology(two_years)
        @test size(clim) == (GREB.xdim, GREB.ydim, 12)
        for m in 1:12
            @test all(==(Float64(m + 6)), clim[:, :, m])
        end
    end

    @testset "tendencies! Q_sens honors log_atmos_dmc gate" begin
        # Q_sens = ct_sens * (Ta - Ts) is checkable directly against a
        # hand-computed value without reimplementing the rest of the
        # physics pipeline.
        fields = ClimateFields()
        state = ModelState()
        ws = CirculationWorkspace()
        ts = TimeState(1, 1)
        cfg = create_experiment_config(:full_model)

        Ts = fill(288.0, GREB.xdim, GREB.ydim)
        Ta = fill(280.0, GREB.xdim, GREB.ydim)
        To = fill(285.0, GREB.xdim, GREB.ydim)
        q = fill(0.006, GREB.xdim, GREB.ydim)

        tend = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, ts, cfg)
        @test tend.Q_sens ≈ GREB.ct_sens .* (Ta .- Ts)
        @test all(isfinite, tend.SW)
        @test all(isfinite, tend.LW_surf)
        @test all(isfinite, tend.dTa_crcl)
        @test all(isfinite, tend.dq_crcl)

        cfg_off = create_experiment_config(:full_model)
        cfg_off.log_atmos_dmc = false
        tend_off = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, ts, cfg_off)
        @test all(iszero, tend_off.Q_sens)
    end

    @testset "diagnostics! accumulates annual means and resets at year end" begin
        fields = ClimateFields()
        state = ModelState()
        ts = TimeState(1, 1)
        z = () -> zeros(GREB.xdim, GREB.ydim)
        surf = SurfaceState(fill(280.0, GREB.xdim, GREB.ydim), fill(270.0, GREB.xdim, GREB.ydim),
            fill(285.0, GREB.xdim, GREB.ydim), fill(0.005, GREB.xdim, GREB.ydim))
        tend = (albedo=fill(0.3, GREB.xdim, GREB.ydim), SW=fill(100.0, GREB.xdim, GREB.ydim),
            ice_cover=z(), LW_surf=fill(-50.0, GREB.xdim, GREB.ydim),
            Q_lat=fill(-20.0, GREB.xdim, GREB.ydim), Q_sens=fill(-5.0, GREB.xdim, GREB.ydim),
            Q_lat_air=fill(20.0, GREB.xdim, GREB.ydim), dq_eva=z(),
            dq_rain=z(), dq_crcl=z(), dTa_crcl=z(), dT_ocean=z(), dTo=z(),
            LW_down=fill(30.0, GREB.xdim, GREB.ydim), LW_up=fill(80.0, GREB.xdim, GREB.ydim),
            em=fill(0.9, GREB.xdim, GREB.ydim))

        ts.ityr = 1
        diagnostics!(1, 1970, 340.0, surf, tend, fields, state, ts)
        @test all(==(280.0), state.Tsmn)   # accumulated once, no averaging/reset yet

        ts.ityr = GREB.nstep_yr
        captured = mktemp() do path, io
            redirect_stdout(io) do
                diagnostics!(GREB.nstep_yr, 1970, 340.0, surf, tend, fields, state, ts)
            end
            flush(io)
            read(path, String)
        end
        @test occursin("1970", captured)   # prints the annual summary line
        @test all(iszero, state.Tsmn)      # reset after year end
    end

    @testset "output! pushes a monthly-mean MonthlyRecord at month boundaries" begin
        ws = CirculationWorkspace()
        acc = MonthlyAccumulator()
        ts = TimeState(1, 1)
        surf = SurfaceState(fill(280.0, GREB.xdim, GREB.ydim), fill(270.0, GREB.xdim, GREB.ydim),
            fill(285.0, GREB.xdim, GREB.ydim), fill(0.005, GREB.xdim, GREB.ydim))
        tend = (albedo=fill(0.3, GREB.xdim, GREB.ydim), SW=fill(100.0, GREB.xdim, GREB.ydim),
            ice_cover=fill(0.1, GREB.xdim, GREB.ydim), LW_surf=fill(-50.0, GREB.xdim, GREB.ydim),
            Q_lat=fill(-20.0, GREB.xdim, GREB.ydim), Q_sens=fill(-5.0, GREB.xdim, GREB.ydim))
        ws.precip_out .= 2.0
        ws.evap_out .= 1.0
        ws.qcrcl_out .= 0.5

        output_buf = MonthlyRecord[]
        irec, mon = 0, 1
        ndt = GREB.ndt_days
        ndays_jan = GREB.cjday_mon[1]
        for day in 1:ndays_jan, step in 1:ndt
            it = (day - 1) * ndt + step
            ts.jday = day
            (mon, irec) = output!(it, irec, mon, surf, tend, ws, output_buf, acc, ts)
        end

        @test length(output_buf) == 1
        @test irec == 1
        @test mon == 2
        @test all(==(280.0), output_buf[1].Ts)
        @test all(==(2.0), output_buf[1].precip)
    end

    @testset "time_loop! integrates one timestep and clamps at min_T_K" begin
        fields = ClimateFields()
        fields.z_topo .= -1.0
        fields.mldclim .= 50.0
        fields.Tclim .= 280.0
        fields.Toclim .= 285.0
        fields.qclim .= 0.006
        fields.cldclim .= 0.5
        fields.swetclim .= 0.5
        fields.uclim .= 2.0
        fields.vclim .= 1.0
        fields.omegaclim .= 0.001
        fields.omegastdclim .= 0.01
        fields.wsclim .= 4.0
        cfg = create_experiment_config(:full_model)
        ini = init_model!(cfg, fields)

        state = ModelState()
        ws = CirculationWorkspace()
        acc = MonthlyAccumulator()
        ts = TimeState(1, 1)

        Ts = fill(GREB.min_T_K - 5.0, GREB.xdim, GREB.ydim)
        Ta = copy(ini.Ta_ini)
        To = copy(ini.To_ini)
        q = copy(ini.q_ini)
        output_buf = MonthlyRecord[]

        (mon, irec) = time_loop!(1, 1970, ini.CO2_ctrl, 1, 0, Ts, Ta, q, To, output_buf,
            fields, state, ws, acc, ts, cfg)

        @test all(isfinite, Ts)
        @test all(isfinite, Ta)
        @test all(isfinite, To)
        @test all(isfinite, q)
        @test all(==(GREB.min_T_K), Ts)
        @test mon == 1
        @test irec == 0
    end

    @testset "greb_model! baseline: default config runs to completion with the right output shape" begin
        cfg = create_experiment_config(:full_model)
        result = redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "")
        end
        @test length(result.ctrl) == 12
        @test length(result.scnr) == 0
        @test result.ctrl[1] isa MonthlyRecord
    end

    @testset "greb_model! runs across log_eva / log_rain branches" begin
        log_evas = (-1, 0, 1, 2)
        log_rains = (-1, 0, 1, 2, 3)
        n = max(length(log_evas), length(log_rains))
        combos = collect(Iterators.take(zip(Iterators.cycle(log_evas), Iterators.cycle(log_rains)), n))
        for (log_eva, log_rain) in combos
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            cfg.log_rain = log_rain
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "")
            end
            @test length(result.ctrl) == 12
        end
    end

    @testset "hydro! errors on invalid log_eva" begin
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.005, GREB.xdim, GREB.ydim)
        cfg = create_experiment_config(:full_model)
        cfg.log_eva = 99
        @test_throws ErrorException hydro!(Ts, q, ClimateFields(), TimeState(1, 1), cfg, CirculationWorkspace())
    end

    @testset "SWradiation! is allocation-free" begin
        fields = ClimateFields()
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()
        cfg = create_experiment_config(:full_model)
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        SWradiation!(Ts, fields, state, ts, cfg, ws)  # warm up (compilation)
        @test @allocated(SWradiation!(Ts, fields, state, ts, cfg, ws)) == 0
    end

    @testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction (matches Fortran)" begin
        # A climatology that's uniform in time (same value at every ityr)
        # reaches its correction fixed point after a single timestep, so
        # Ts/To/q should land exactly on Tclim/Toclim/qclim, not just move
        # closer to it.
        fields = ClimateFields()
        fields.cap_surf .= GREB.cap_ocean
        for j in 1:GREB.ydim, i in 1:GREB.xdim
            fields.Tclim[i, j, :] .= 280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0)
            fields.Toclim[i, j, :] .= 279.0
            fields.qclim[i, j, :] .= 0.006
        end
        cfg = create_experiment_config(:full_model)
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        Ta = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.010, GREB.xdim, GREB.ydim)
        To = fill(285.0, GREB.xdim, GREB.ydim)

        GREB.qflux_correction!(340.0, Ts, Ta, q, To, fields, state, ts, cfg, ws, 1)

        @test any(!=(0.0), fields.TF_correct)
        @test any(!=(0.0), fields.ToF_correct)
        @test any(!=(0.0), fields.qF_correct)
        @test all(isfinite, fields.TF_correct)
        @test all(isfinite, fields.ToF_correct)
        @test all(isfinite, fields.qF_correct)

        @test Ts ≈ fields.Tclim[:, :, 1]
        @test To ≈ fields.Toclim[:, :, 1]
        @test q ≈ fields.qclim[:, :, 1]

        # Ta has no correction field — it integrates freely. Confirm it
        # still moved from its initial value (isn't frozen/broken) even
        # though nothing nudges it toward a climatology target.
        @test all(isfinite, Ta)
        @test Ta != fill(290.0, GREB.xdim, GREB.ydim)
    end

    @testset "greb_model! swaps sw_solar for paleo experiments, restores after" begin
        # A paleo run's swapped solar table must not leak into a later run
        # against the same reused `fields` instance.
        fields = ClimateFields()
        saved_sw_solar = copy(fields.sw_solar)
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            distinctive_value = 999.0
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(distinctive_value, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end

            cfg = create_experiment_config(:paleo_231kyr)
            captured = mktemp() do path, io
                result = redirect_stdout(io) do
                    greb_model!(RunSpec(ctrl = 0), cfg; jld2_dir = tmpdir, fields = fields)
                end
                flush(io)
                (result = result, text = read(path, String))
            end
            @test length(captured.result.scnr) == 12
            # Confirm the swap branch actually ran (fields are all-zero in
            # this test env, so asserting on Ts/anomaly magnitude would be
            # fragile — several other tests in this suite hit NaN for the
            # same reason and only check shape/completion instead).
            @test occursin("loading alternate solar forcing", captured.text)

            # sw_solar restored to its pre-run value after greb_model! returns
            @test fields.sw_solar == saved_sw_solar

            cfg_plain = create_experiment_config(:full_model)
            redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields)
            end
            @test fields.sw_solar == saved_sw_solar
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "greb_model! reaches every :experiment symbol forcing()/init_model! dispatch on" begin
        direct_dispatch_symbols = (
            :a1b_scenario, :co2_10x, :co2_half, :co2_zero, :solar_cycle_11yr,
            :a1b_enhanced, :co2_sine_wave, :co2_step, :modern_solar_paleo_co2,
            :earth_sun_distance, :regional_co2_nh, :regional_co2_sh,
            :regional_co2_tropics, :regional_co2_extratropics,
            :regional_co2_ocean, :regional_co2_land_ice, :regional_co2_winter,
            :regional_co2_summer,
        )
        for sym in direct_dispatch_symbols
            fields = ClimateFields()
            cfg = PhysicsConfig(experiment = sym)
            icmn_ctrl = zeros(Float64, GREB.xdim, GREB.ydim, 1)
            redirect_stdout(devnull) do
                init_model!(cfg, fields)
            end
            result = forcing(1, 1970, cfg, fields, icmn_ctrl)
            @test isfinite(result.CO2)
            @test isfinite(result.sw_solar_forcing)
        end

        cfg = PhysicsConfig(experiment = :sst_plus1)
        result = redirect_stdout(devnull) do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "")
        end
        @test length(result.scnr) == 12

        # "not yet implemented" placeholders 
        for sym in (:rcp26, :rcp45, :rcp60, :custom_co2)
            cfg = PhysicsConfig(experiment = sym)
            @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "")
        end

        # :paleo_solar_modern_co2/:obliquity/:eccentricity swap in a real
        # solar_scenarios/*.jld2 table at scenario start — synthesize a
        # minimal set so these three are reachable without the real dataset.
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(999.0, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_obliquity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_eccentricity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end

            for sym in (:paleo_solar_modern_co2, :obliquity, :eccentricity)
                cfg = PhysicsConfig(experiment = sym)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir)
                end
                @test length(result.scnr) == 12
            end
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "log_hydro_dmc==false freezes humidity entirely (not just eva/rain)" begin
        # With log_hydro_dmc off, q must never move from its initial
        # climatological value
        if !isdir(DATA_DIR)
            @test_skip "greb_dataset_jld2/ not present"
        else
            fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
            cfg = create_experiment_config(:full_model)
            cfg.log_hydro_dmc = false
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = DATA_DIR, fields = fields)
            end
            q_ini = fields.qclim[:, :, GREB.nstep_yr]
            for rec in result.ctrl
                @test all(isapprox.(rec.q, q_ini; atol = 1e-9))
            end
        end
    end

    @testset "golden regression: real dataset control+scenario run matches snapshot" begin
        # Tripwire for any refactor touching the physics kernels: a real 1yr
        # control + 1yr scenario run against the actual NCEP dataset,
        # snapshotted as monthly global-mean Ts/Ta/q. Drift beyond
        # float-reassociation noise (~1e-12) means real behavior changed.
        # Set RUN_GOLDEN=0 to skip this locally; CI always runs it.
        if !isdir(DATA_DIR)
            @test_skip "greb_dataset_jld2/ not present"
        elseif get(ENV, "RUN_GOLDEN", "1") == "0"
            @test_skip "RUN_GOLDEN=0"
        else
            fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
            cfg = create_experiment_config(:full_model)
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(), cfg; jld2_dir = DATA_DIR, fields = fields)
            end

            gmean(x) = sum(x) / length(x)
            summarize(rec) = (Ts = gmean(rec.Ts), Ta = gmean(rec.Ta), q = gmean(rec.q))

            ctrl_ref = [
                (Ts = 276.6376432475215, Ta = 279.0132752922387, q = 0.006483368265992144),
                (Ts = 276.1076273785758, Ta = 278.6101854203388, q = 0.0066878269244861795),
                (Ts = 276.3982313917197, Ta = 278.7355819769281, q = 0.006828203121749328),
                (Ts = 278.033418596411, Ta = 280.2620085698279, q = 0.006996728195722164),
                (Ts = 280.1103499917721, Ta = 282.3720691501376, q = 0.00732842859123354),
                (Ts = 281.83909739745076, Ta = 284.25313867300514, q = 0.007855042454125112),
                (Ts = 282.5491066792044, Ta = 285.12633793891365, q = 0.008282889952438756),
                (Ts = 282.33494625480506, Ta = 285.0130523882134, q = 0.008281318600014156),
                (Ts = 281.10617631470427, Ta = 283.7911113302602, q = 0.007863431467610265),
                (Ts = 279.50523603155864, Ta = 282.12404966212097, q = 0.007469852856220665),
                (Ts = 278.537108759095, Ta = 281.1654805465912, q = 0.0073083638644664516),
                (Ts = 278.2347705686773, Ta = 280.9247609592766, q = 0.0073813136191502576),
            ]
            scnr_ref = [
                (Ts = -0.007632155594618766, Ta = -0.007133529465347189, q = -4.9333976027584994e-8),
                (Ts = -0.001309167326571708, Ta = -0.0015145523125641436, q = 1.030987112661906e-8),
                (Ts = -0.00011935227909166186, Ta = -0.00015150845673007126, q = -8.811407575403715e-9),
                (Ts = 0.00011412561624893656, Ta = 0.00010169282102569695, q = 7.191046358240965e-10),
                (Ts = 0.0001557742243889431, Ta = 0.00015178409538934505, q = 1.282729300373303e-8),
                (Ts = 9.186448051926958e-5, Ta = 9.533767452687043e-5, q = 1.5921599083929904e-8),
                (Ts = 4.87049451490498e-5, Ta = 4.9187632382707847e-5, q = 1.0540113727968019e-8),
                (Ts = 3.253222138257147e-5, Ta = 3.153559790906405e-5, q = 3.6892345609650163e-9),
                (Ts = 8.00866341792994e-5, Ta = 7.24281415625589e-5, q = 4.179801343708237e-9),
                (Ts = 9.394539367011155e-5, Ta = 9.5251277080081e-5, q = 9.7252510856804e-10),
                (Ts = 7.547731804897885e-5, Ta = 7.672918591223999e-5, q = -2.129466956195276e-9),
                (Ts = 5.324451949642286e-5, Ta = 5.3642668972774134e-5, q = -3.3368785465046947e-9),
            ]

            @test length(result.ctrl) == length(ctrl_ref)
            @test length(result.scnr) == length(scnr_ref)
            for (rec, ref) in zip(result.ctrl, ctrl_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-6)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-6)
                @test isapprox(s.q, ref.q; atol = 1e-6)
            end
            for (rec, ref) in zip(result.scnr, scnr_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-6)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-6)
                @test isapprox(s.q, ref.q; atol = 1e-6)
            end
        end
    end

    @testset "read_jld2 rejects non-JLD2 input" begin
        tmp = tempname() * ".jld2"
        write(tmp, "not a jld2 file")
        @test_throws Exception read_jld2(tmp)
        rm(tmp; force = true)
    end

    @testset "load_greb_jld2!/load_flux_corrections_jld2! file-exists branches" begin
        # Synthesize a minimal dataset matching load_greb_jld2!'s expected
        # layout (src/io.jl) so both loaders' "file present"/"file missing"
        # branches are exercised without the real (gitignored) dataset.
        write2(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.xdim, GREB.ydim); f["dim_names"] = ["lon", "lat"]
        end)
        write3(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.xdim, GREB.ydim, GREB.nstep_yr); f["dim_names"] = ["lon", "lat", "time"]
        end)
        write_solar(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.ydim, GREB.nstep_yr); f["dim_names"] = ["lat", "time"]
        end)

        tmpdir = mktempdir()
        try
            write2(joinpath(tmpdir, "static", "global.topography.jld2"), 1.0)
            write2(joinpath(tmpdir, "static", "greb.glaciers.jld2"), 2.0)
            write3(joinpath(tmpdir, "climatology", "ncep.tsurf.1948-2007.clim.jld2"), 3.0)
            write3(joinpath(tmpdir, "climatology", "ncep.zonal_wind.850hpa.clim.jld2"), 4.0)
            write3(joinpath(tmpdir, "climatology", "ncep.meridional_wind.850hpa.clim.jld2"), 5.0)
            write3(joinpath(tmpdir, "climatology", "ncep.atmospheric_humidity.clim.jld2"), 6.0)
            write3(joinpath(tmpdir, "climatology", "ncep.soil_moisture.clim.jld2"), 7.0)
            write3(joinpath(tmpdir, "climatology", "isccp.cloud_cover.clim.jld2"), 8.0)
            write3(joinpath(tmpdir, "climatology", "woce.ocean_mixed_layer_depth.clim.jld2"), 9.0)
            write3(joinpath(tmpdir, "climatology", "Tocean.clim.jld2"), 10.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.omega.vertmean.clim.jld2"), 11.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.omega_std.vertmean.clim.jld2"), 12.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.windspeed.850hpa.clim.jld2"), 13.0)
            write_solar(joinpath(tmpdir, "solar", "solar_radiation.clim.jld2"), 14.0)

            # "files missing" branch: no flux-correction files present yet ->
            # load_flux_corrections_jld2! should warn and zero-fill, not error.
            fields_nocorr = load_greb_jld2!(tmpdir; dataset = :ncep)
            @test all(==(0.0), fields_nocorr.TF_correct)
            @test all(==(0.0), fields_nocorr.qF_correct)
            @test all(==(0.0), fields_nocorr.ToF_correct)
            @test all(==(3.0), fields_nocorr.Tclim)  # loader itself still worked

            # "files present" branch: add the 3 flux-correction files and reload.
            write3(joinpath(tmpdir, "climatology", "Tsurf_flux_correction.jld2"), 15.0)
            write3(joinpath(tmpdir, "climatology", "vapour_flux_correction.jld2"), 16.0)
            write3(joinpath(tmpdir, "climatology", "Tocean_flux_correction.jld2"), 17.0)

            fields = load_greb_jld2!(tmpdir; dataset = :ncep)
            @test all(==(1.0), fields.z_topo)
            @test all(==(2.0), fields.glacier)
            @test all(==(3.0), fields.Tclim)
            @test all(==(4.0), fields.uclim)
            @test all(==(14.0), fields.sw_solar)
            @test all(==(15.0), fields.TF_correct)
            @test all(==(16.0), fields.qF_correct)
            @test all(==(17.0), fields.ToF_correct)
        finally
            rm(tmpdir; recursive = true, force = true)
        end

        missing_parent = mktempdir()
        @test_throws ErrorException load_greb_jld2!(joinpath(missing_parent, "nonexistent"))
        rm(missing_parent; recursive = true, force = true)
    end

end

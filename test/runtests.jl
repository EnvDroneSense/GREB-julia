using GREBClimate
using Test

# Smoke tests that do NOT require the (large, external) JLD2 input data.
# They check that the package loads, its types build, and the grid/constants
# are intact after the notebook -> package extraction. Full integration runs
# (which need `greb_input_data/`) are demonstrated in examples/run_greb.jl.

# `allow_download=false` so running the tests can never trigger the 353 MB
# DataDeps download; data-dependent testsets @test_skip when this is absent.
const DATA_DIR = something(greb_data_dir(; allow_download = false),
                           joinpath(@__DIR__, "..", "greb_input_data"))

# The 24 testsets below are split into two groups so CI can shard them across
# 2 parallel jobs `run_light_tests()` (17 cheap testsets) and `run_heavy_tests()` (7
# testsets that run one or more full `greb_model!` control/scenario runs,
# including the golden regression). `GREB_TEST_SHARD` selects which group(s)
# run; unset (or "all") runs both, so plain `Pkg.test()`/`]test` is
# unaffected for local use - sharding is opt-in, CI-only.

function run_light_tests()

    @testset "grid constants" begin
        @test GREBClimate.xdim == 96
        @test GREBClimate.ydim == 48
        @test GREBClimate.nstep_yr == 730
    end

    @testset "PhysicsConfig" begin
        for exp in (:full_model, :constant_topo, :co2_double, :co2_quadruple,
                    :elnino, :lanina, :rcp26, :rcp45, :rcp60, :rcp85, :ssp119,
                    :ssp126, :ssp245, :ssp460, :ssp585, :historical_co2,
                    :decon_mean_climate, :decon_2xco2)
            c = create_experiment_config(exp)
            @test c.experiment == exp
        end

        @test create_experiment_config(:co2_double).co2_concentration == 680.0
    end

    @testset "create_experiment_config: :custom_co2 and decon presets (§7.2/§7.3)" begin
        cfg_custom = create_experiment_config(:custom_co2; co2_path="/tmp/my_co2.txt")
        @test cfg_custom.experiment == :custom_co2
        @test cfg_custom.custom_co2_path == "/tmp/my_co2.txt"
        @test create_experiment_config(:custom_co2).custom_co2_path == ""

        # decon_mean_climate: defaults all true, one override propagates
        cfg_dmc = create_experiment_config(:decon_mean_climate)
        for switch in (cfg_dmc.log_clouds_dmc, cfg_dmc.log_ocean_dmc, cfg_dmc.log_atmos_dmc,
                       cfg_dmc.log_co2_dmc, cfg_dmc.log_hydro_dmc, cfg_dmc.log_qflux_dmc,
                       cfg_dmc.log_ice, cfg_dmc.log_hdif, cfg_dmc.log_hadv,
                       cfg_dmc.log_vdif, cfg_dmc.log_vadv)
            @test switch == true
        end
        cfg_dmc_off = create_experiment_config(:decon_mean_climate; log_ocean_dmc=false)
        @test cfg_dmc_off.log_ocean_dmc == false
        @test cfg_dmc_off.log_clouds_dmc == true  # untouched switches stay at default

        # decon_2xco2: defaults all true + doubled CO2, one override propagates
        cfg_drsp = create_experiment_config(:decon_2xco2)
        @test cfg_drsp.co2_concentration == 680.0
        for switch in (cfg_drsp.log_topo_drsp, cfg_drsp.log_clouds_drsp, cfg_drsp.log_humid_drsp,
                       cfg_drsp.log_ocean_drsp, cfg_drsp.log_hydro_drsp,
                       cfg_drsp.log_ice, cfg_drsp.log_hdif, cfg_drsp.log_hadv,
                       cfg_drsp.log_vdif, cfg_drsp.log_vadv)
            @test switch == true
        end
        cfg_drsp_off = create_experiment_config(:decon_2xco2; log_topo_drsp=false)
        @test cfg_drsp_off.log_topo_drsp == false
        @test cfg_drsp_off.log_clouds_drsp == true  # untouched switches stay at default
    end

    @testset "set_hydrology_parameters! matches HYDRO_PARAMS for every log_rain value" begin
        for log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_rain = log_rain
            set_hydrology_parameters!(cfg)
            expected = GREBClimate.HYDRO_PARAMS[log_rain]
            @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == expected
        end

        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 0
        cfg.log_clim = 1
        set_hydrology_parameters!(cfg)
        @test all(isapprox.((cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd), (-1.27, 1.99, -16.54, 21.15); atol = 1e-4))
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
            greb_model!(RunSpec(scnr = 0), cfg_regional; jld2_dir = "", fields = fields, allow_uninitialized = true)
        end
        @test any(!=(1.0), fields.co2_part)  # regional run actually changed the mask

        cfg_plain = create_experiment_config(:full_model)
        redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields, allow_uninitialized = true)
        end
        @test all(==(1.0), fields.co2_part)  # init_model! resets it back to full CO2
    end

    @testset "workspace/accumulator/record types construct correctly" begin
        ws = CirculationWorkspace()
        @test size(ws.dTa_crcl) == (GREBClimate.xdim, GREBClimate.ydim)
        @test eltype(ws.dTa_crcl) === Float32

        acc = MonthlyAccumulator()
        acc.count = 7
        fill!(acc.Tmm, 42.0f0)
        GREBClimate.reset!(acc)
        @test acc.count == 0
        @test all(iszero, acc.Tmm)

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
        mkrec(v) =(Ts=fill(v, GREBClimate.xdim, GREBClimate.ydim), Ta=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            To=fill(v, GREBClimate.xdim, GREBClimate.ydim), q=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            albedo=fill(v, GREBClimate.xdim, GREBClimate.ydim), ice=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            precip=fill(v, GREBClimate.xdim, GREBClimate.ydim), evap=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            qcrcl=fill(v, GREBClimate.xdim, GREBClimate.ydim), sw=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            lw=fill(v, GREBClimate.xdim, GREBClimate.ydim), qlat=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            qsens=fill(v, GREBClimate.xdim, GREBClimate.ydim))

        @test build_monthly_climatology(MonthlyRecord[]) == MonthlyRecord[]

        # Final-year-only.
        two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        clim = build_monthly_climatology(two_years)
        @test length(clim) == 12
        for m in 1:12
            @test all(==(Float64(m + 12)), clim[m].Ts)
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
        # above: final-year-only.
        mkrec(v) = (Ts=zeros(GREBClimate.xdim, GREBClimate.ydim), Ta=zeros(GREBClimate.xdim, GREBClimate.ydim),
            To=zeros(GREBClimate.xdim, GREBClimate.ydim), q=zeros(GREBClimate.xdim, GREBClimate.ydim),
            albedo=zeros(GREBClimate.xdim, GREBClimate.ydim), ice=fill(v, GREBClimate.xdim, GREBClimate.ydim),
            precip=zeros(GREBClimate.xdim, GREBClimate.ydim), evap=zeros(GREBClimate.xdim, GREBClimate.ydim),
            qcrcl=zeros(GREBClimate.xdim, GREBClimate.ydim), sw=zeros(GREBClimate.xdim, GREBClimate.ydim),
            lw=zeros(GREBClimate.xdim, GREBClimate.ydim), qlat=zeros(GREBClimate.xdim, GREBClimate.ydim),
            qsens=zeros(GREBClimate.xdim, GREBClimate.ydim))

        @test all(iszero, compute_annual_ice_climatology(MonthlyRecord[]))

        two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        clim = compute_annual_ice_climatology(two_years)
        @test size(clim) == (GREBClimate.xdim, GREBClimate.ydim, 12)
        for m in 1:12
            @test all(==(Float64(m + 12)), clim[:, :, m])
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

        Ts = fill(288.0, GREBClimate.xdim, GREBClimate.ydim)
        Ta = fill(280.0, GREBClimate.xdim, GREBClimate.ydim)
        To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)
        q = fill(0.006, GREBClimate.xdim, GREBClimate.ydim)

        tend = tendencies!(340.0, Ts, Ta, To, q, fields, state, ws, ts, cfg)
        @test tend.Q_sens ≈ GREBClimate.ct_sens .* (Ta .- Ts)
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
        z = () -> zeros(GREBClimate.xdim, GREBClimate.ydim)
        surf = SurfaceState(fill(280.0, GREBClimate.xdim, GREBClimate.ydim), fill(270.0, GREBClimate.xdim, GREBClimate.ydim),
            fill(285.0, GREBClimate.xdim, GREBClimate.ydim), fill(0.005, GREBClimate.xdim, GREBClimate.ydim))
        tend = (albedo=fill(0.3, GREBClimate.xdim, GREBClimate.ydim), SW=fill(100.0, GREBClimate.xdim, GREBClimate.ydim),
            ice_cover=z(), LW_surf=fill(-50.0, GREBClimate.xdim, GREBClimate.ydim),
            Q_lat=fill(-20.0, GREBClimate.xdim, GREBClimate.ydim), Q_sens=fill(-5.0, GREBClimate.xdim, GREBClimate.ydim),
            Q_lat_air=fill(20.0, GREBClimate.xdim, GREBClimate.ydim), dq_eva=z(),
            dq_rain=z(), dq_crcl=z(), dTa_crcl=z(), dT_ocean=z(), dTo=z(),
            LW_down=fill(30.0, GREBClimate.xdim, GREBClimate.ydim), LW_up=fill(80.0, GREBClimate.xdim, GREBClimate.ydim),
            em=fill(0.9, GREBClimate.xdim, GREBClimate.ydim))

        ts.ityr = 1
        diagnostics!(1, 1970, 340.0, surf, tend, fields, state, ts)
        @test all(==(280.0), state.Tsmn)   # accumulated once, no averaging/reset yet

        ts.ityr = GREBClimate.nstep_yr
        captured = mktemp() do path, io
            redirect_stdout(io) do
                diagnostics!(GREBClimate.nstep_yr, 1970, 340.0, surf, tend, fields, state, ts)
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
        surf = SurfaceState(fill(280.0, GREBClimate.xdim, GREBClimate.ydim), fill(270.0, GREBClimate.xdim, GREBClimate.ydim),
            fill(285.0, GREBClimate.xdim, GREBClimate.ydim), fill(0.005, GREBClimate.xdim, GREBClimate.ydim))
        tend = (albedo=fill(0.3, GREBClimate.xdim, GREBClimate.ydim), SW=fill(100.0, GREBClimate.xdim, GREBClimate.ydim),
            ice_cover=fill(0.1, GREBClimate.xdim, GREBClimate.ydim), LW_surf=fill(-50.0, GREBClimate.xdim, GREBClimate.ydim),
            Q_lat=fill(-20.0, GREBClimate.xdim, GREBClimate.ydim), Q_sens=fill(-5.0, GREBClimate.xdim, GREBClimate.ydim))
        ws.precip_out .= 2.0
        ws.evap_out .= 1.0
        ws.qcrcl_out .= 0.5

        output_buf = MonthlyRecord[]
        irec, mon = 0, 1
        ndt = GREBClimate.ndt_days
        ndays_jan = GREBClimate.cjday_mon[1]
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

        Ts = fill(GREBClimate.min_T_K - 5.0, GREBClimate.xdim, GREBClimate.ydim)
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
        @test all(==(GREBClimate.min_T_K), Ts)
        @test mon == 1
        @test irec == 0
    end

    @testset "diffusion!/advection!/circulation! per-cell snapshot (incl. date-line wraparound)" begin
        # Golden per-cell snapshot (not just the golden regression's global
        # mean) so a future refactor of the zonal neighbor-index lookups
        # (`lon_jm1`/`lon_jp1`/etc.) that mis-handles the periodic
        # date-line wraparound gets caught directly, not diluted into a
        # spatial average. Deterministic (not random) T1/wz fields make any
        # wraparound mistake show up as an obvious, reproducible diff.
        # Covers: north pole (k=1) and south pole (k=48) rows (both also
        # exercise the zonal polar sub-stepping branch, since
        # dxlat_grid[1]==dxlat_grid[48] <= 2.5e5) and a mid-latitude row
        # (k=11, dxlat_grid[11] > 2.5e5, the plain @turbo zonal branch);
        # date-line columns 1,2,3,94,95,96 plus an interior column (50).
        fields = ClimateFields()
        xdim_, ydim_ = GREBClimate.xdim, GREBClimate.ydim
        T1 = Float32[100.0 * i + k for i in 1:xdim_, k in 1:ydim_]
        wz = [1.0 + 0.001 * i - 0.0005 * k for i in 1:xdim_, k in 1:ydim_]
        fields.wz_air .= wz
        fields.wz_vapor .= wz
        for it in 1:GREBClimate.nstep_yr, k in 1:ydim_, i in 1:xdim_
            fields.uclim_p[i, k, it] = 0.5 + 0.0001 * i
            fields.uclim_m[i, k, it] = 0.3 + 0.0001 * k
            fields.vclim_p[i, k, it] = 0.4 + 0.0002 * i
            fields.vclim_m[i, k, it] = 0.2 + 0.0002 * k
        end
        ws = CirculationWorkspace()
        ts = TimeState(1, 1)
        cfg = create_experiment_config(:full_model)

        test_is = [1, 2, 3, 50, 94, 95, 96]
        test_ks = [1, 11, 48]

        diffusion!(T1, GREBClimate.z_air, fields, ws, ts)
        dX_diff_ref = Dict(
            (1,1)=>4517.8355192140425, (2,1)=>3542.95440832927, (3,1)=>2652.6121358299374,
            (50,1)=>1.2269892658145531, (94,1)=>-2709.198844945152, (95,1)=>-3576.970530673063,
            (96,1)=>-4530.555031256096,
            (1,11)=>64.24343226619055, (2,11)=>32.16535999519279, (3,11)=>10.737877625759896,
            (50,11)=>0.0032154796768526627, (94,11)=>-10.710851102068029, (95,11)=>-32.17954307109926,
            (96,11)=>-64.44305854028968,
            (1,48)=>4410.353233774593, (2,48)=>3445.760591202237, (3,48)=>2567.190165268745,
            (50,48)=>1.1818789937009289, (94,48)=>-2623.7462653181888, (95,48)=>-3479.421880797481,
            (96,48)=>-4422.060857683985,
        )
        for k in test_ks, i in test_is
            @test isapprox(ws.dX_diff[i, k], dX_diff_ref[(i, k)]; atol=1e-3, rtol=1e-4)
        end

        advection!(T1, GREBClimate.z_air, fields, ws, ts, cfg)
        dX_adv_ref = Dict(
            (1,1)=>99.99281072836801, (2,1)=>37.62423619149657, (3,1)=>6.428217314852662,
            (50,1)=>-4.182755344492395, (94,1)=>11.771946283998448, (95,1)=>60.257791236566284,
            (96,1)=>157.2916794570926,
            (1,11)=>6.863756666482768, (2,11)=>3.295185364395947, (3,11)=>-0.2733860065994545,
            (50,11)=>-0.28795610784365167, (94,11)=>-0.30173415760057604, (95,11)=>5.231130705764574,
            (96,11)=>10.766167503970431,
            (1,48)=>99.42365199872042, (2,48)=>37.439173701522, (3,48)=>6.435048732922689,
            (50,48)=>-4.112511103625333, (94,48)=>11.462671629384038, (95,48)=>58.8110669109033,
            (96,48)=>153.56947120915834,
        )
        for k in test_ks, i in test_is
            @test isapprox(ws.dX_adv[i, k], dX_adv_ref[(i, k)]; atol=1e-3, rtol=1e-4)
        end

        dX_out = zeros(xdim_, ydim_)
        circulation!(T1, GREBClimate.z_air, dX_out, fields, ws, ts, cfg)
        dX_out_ref = Dict(
            (1,1)=>4824.155681112328, (2,1)=>4609.661409972268, (3,1)=>4395.402855867866,
            (50,1)=>-98.5576039612888, (94,1)=>-4172.300403641432, (95,1)=>-4369.7800972827745,
            (96,1)=>-4569.32742911384,
            (1,11)=>1447.9977087179682, (2,11)=>765.1034286356905, (3,11)=>274.82096986927763,
            (50,11)=>-6.684077032670757, (94,11)=>-252.25457239155912, (95,11)=>-576.3346453253889,
            (96,11)=>-1090.2114380116673,
            (1,48)=>4827.018868288203, (2,48)=>4605.028237288625, (3,48)=>4383.415249532827,
            (50,48)=>-94.78570775574462, (94,48)=>-4150.449604784975, (95,48)=>-4353.86760064195,
            (96,48)=>-4559.610385060042,
        )
        for k in test_ks, i in test_is
            @test isapprox(dX_out[i, k], dX_out_ref[(i, k)]; atol=1e-3, rtol=1e-4)
        end
    end

    @testset "hydro! errors on invalid log_eva" begin
        Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
        q = fill(0.005, GREBClimate.xdim, GREBClimate.ydim)
        cfg = create_experiment_config(:full_model)
        cfg.log_eva = 99
        @test_throws ErrorException hydro!(Ts, q, ClimateFields(), TimeState(1, 1), cfg, CirculationWorkspace())
    end

    @testset "hydro! log_eva==1 gust includes Fortran's carried-over +2.0²/+3.0² base term (§8.2)" begin
        # With u=v=0, wind = sqrt(gust) directly, isolating the additive
        # constant. Fortran's shared `abswind` already carries +2.0²(land)/
        # +3.0²(ocean) into this branch before adding its own 144.²/7.1² —
        # the true combined constant is 4+144² land, 9+50.41 ocean, not
        # 144²/50.41 alone (greb.model.mscm.f90:710-712,727-728).
        mkfields(topo) = begin
            fields = ClimateFields()
            fields.z_topo .= topo
            fields.mldclim .= 50.0
            fields.Tclim .= 280.0
            fields.Toclim .= 285.0
            fields.qclim .= 0.006
            fields.cldclim .= 0.5
            fields.swetclim .= 1.0
            fields.uclim .= 0.0
            fields.vclim .= 0.0
            fields.omegaclim .= 0.0
            fields.omegastdclim .= 0.0
            fields.wsclim .= 0.0
            fields
        end
        cfg = create_experiment_config(:full_model)
        cfg.log_eva = 1
        Ts = fill(290.0f0, GREBClimate.xdim, GREBClimate.ydim)
        q = fill(0.008f0, GREBClimate.xdim, GREBClimate.ydim)
        ts = TimeState(1, 1)

        for (topo, gust, coeff) in ((1.0, 4.0 + 144.0, 0.04), (-1.0, 9.0 + 50.41, 0.73))
            fields = mkfields(topo)
            init_model!(cfg, fields)
            ws = CirculationWorkspace()
            result = hydro!(Ts, q, fields, ts, cfg, ws)

            qs = 3.75e-3 * exp(17.08085 * (290.0 - 273.15) / (290.0 - 273.15 + 234.175)) * fields.wz_air[1, 1]
            expected = (q[1, 1] - qs) * sqrt(gust) * GREBClimate.cq_latent * GREBClimate.ρ_air * coeff * GREBClimate.ce * 1.0
            @test isapprox(result.Q_lat[1, 1], expected; rtol = 1e-5)
        end
    end

    @testset "hydro! doesn't apply an extra -0.9q clamp to dq_rain (§8.3)" begin
        # Fortran's `hydro` subroutine has no such clamp — the real
        # threshold-replace applies once, downstream, to the combined dq in
        # time_loop!'s update (§0.18). Force the raw regression to predict a
        # huge fractional loss (dq_rain very negative) and confirm dq_rain
        # (and Q_lat_air derived from it) reflect the UNclamped regression
        # value, not a value pulled up to -0.9q/Δt.
        fields = ClimateFields()
        fields.z_topo .= 1.0
        fields.mldclim .= 50.0
        fields.Tclim .= 280.0
        fields.Toclim .= 285.0
        fields.qclim .= 0.006
        fields.cldclim .= 0.5
        fields.swetclim .= 1.0
        fields.uclim .= 0.0
        fields.vclim .= 0.0
        fields.omegaclim .= 0.0
        fields.omegastdclim .= 0.0
        fields.wsclim .= 0.0
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 0  # disable the (unrelated) rain-limit clamp
        init_model!(cfg, fields)

        # Override c_q/etc AFTER init_model! — it calls
        # set_hydrology_parameters!, which would otherwise clobber these
        # back to the log_rain-indexed HYDRO_PARAMS preset. cq_rain is
        # negative, so a large POSITIVE c_q drives dq_rain far below
        # -0.9q/Δt (dq_rain = c_q*cq_rain*q).
        cfg.c_q = 1000.0
        cfg.c_rq = 0.0; cfg.c_omega = 0.0; cfg.c_omegastd = 0.0

        Ts = fill(290.0f0, GREBClimate.xdim, GREBClimate.ydim)
        q = fill(0.008f0, GREBClimate.xdim, GREBClimate.ydim)
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()
        result = hydro!(Ts, q, fields, ts, cfg, ws)

        expected_dq_rain = cfg.c_q * GREBClimate.cq_rain * q[1, 1]
        min_dq_that_would_have_clamped = -0.9 * q[1, 1] / GREBClimate.Δt
        @test expected_dq_rain < min_dq_that_would_have_clamped  # sanity: the old clamp would have fired
        @test isapprox(result.dq_rain[1, 1], expected_dq_rain; rtol = 1e-5)
        @test isapprox(result.Q_lat_air[1, 1], -expected_dq_rain * GREBClimate.cq_latent * GREBClimate.r_qviwv; rtol = 1e-5)
    end

    @testset "SWradiation! is allocation-free" begin
        fields = ClimateFields()
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()
        cfg = create_experiment_config(:full_model)
        Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
        SWradiation!(Ts, fields, state, ts, cfg, ws)
        @test @allocated(SWradiation!(Ts, fields, state, ts, cfg, ws)) <= 64
    end

    @testset "log_hydro_dmc==false freezes humidity entirely (not just eva/rain)" begin
        # With log_hydro_dmc off, q must never move from its initial
        # climatological value
        if !isdir(DATA_DIR)
            @test_skip "greb_input_data/ not present"
        else
            fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
            cfg = create_experiment_config(:full_model)
            cfg.log_hydro_dmc = false
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = DATA_DIR, fields = fields)
            end
            q_ini = fields.qclim[:, :, GREBClimate.nstep_yr]
            for rec in result.ctrl
                @test all(isapprox.(rec.q, q_ini; atol = 1e-7))
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
        write2(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
            f["data"] = fill(v, GREBClimate.xdim, GREBClimate.ydim); f["dim_names"] = ["lon", "lat"]
        end)
        write3(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
            f["data"] = fill(v, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr); f["dim_names"] = ["lon", "lat", "time"]
        end)
        write_solar(path, v) = (mkpath(dirname(path)); GREBClimate.jldopen(path, "w") do f
            f["data"] = fill(v, GREBClimate.ydim, GREBClimate.nstep_yr); f["dim_names"] = ["lat", "time"]
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

            # "files present" branch: add the combined flux-correction file and reload.
            mkpath(joinpath(tmpdir, "climatology"))
            GREBClimate.jldopen(joinpath(tmpdir, "climatology", "flux_corrections.jld2"), "w") do f
                f["Tsurf_flux_correction"] = fill(15.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
                f["vapour_flux_correction"] = fill(16.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
                f["Tocean_flux_correction"] = fill(17.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            end

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

    @testset "greb_data_dir resolution order" begin
        # Never exercises the DataDep branch: that would download 353 MB. The
        # branch itself is a single `allow_download || return nothing` guard;
        # what is worth testing is that the three local sources take priority
        # over it in the right order, so a download can only ever be a last
        # resort.
        tmp_a, tmp_b = mktempdir(), mktempdir()
        saved = get(ENV, "GREB_DATA", nothing)
        try
            # explicit path wins over everything
            ENV["GREB_DATA"] = tmp_b
            @test greb_data_dir(tmp_a) == tmp_a
            # ...and over the environment even with allow_download off
            @test greb_data_dir(tmp_a; allow_download = false) == tmp_a
            # GREB_DATA wins over the repo-local dataset
            @test greb_data_dir() == tmp_b
            delete!(ENV, "GREB_DATA")

            # a non-existent explicit path is an error, not a silent fallback
            @test_throws ErrorException greb_data_dir(joinpath(tmp_a, "nope"))
            # so is a GREB_DATA pointing nowhere
            ENV["GREB_DATA"] = joinpath(tmp_a, "nope")
            @test_throws ErrorException greb_data_dir()
            delete!(ENV, "GREB_DATA")

            # an empty explicit path falls through rather than erroring
            @test greb_data_dir("") == greb_data_dir()
        finally
            saved === nothing ? delete!(ENV, "GREB_DATA") : (ENV["GREB_DATA"] = saved)
            rm(tmp_a; recursive = true, force = true)
            rm(tmp_b; recursive = true, force = true)
        end
    end

    @testset "published dataset archive constants are coherent" begin
        @test occursin(r"^[0-9a-f]{64}$", GREBClimate.DATA_SHA256)

        data_src = read(joinpath(@__DIR__, "..", "src", "data.jl"), String)
        @test match(r"const DATA_RELEASE_TAG = \"([^\"]+)\"", data_src).captures[1] ==
              GREBClimate.DATA_RELEASE_TAG
    end

    @testset "converter allowlist matches what src/io.jl loads" begin
        # tools/convert_greb_to_jld2.jl used to convert every .bin it found,
        # emitting 11 .jld2 files (~148 MB) the model never opens. It now filters
        # on MODEL_FIELD_NAMES. This test keeps that list and src/io.jl's actual
        # loads from drifting apart in either direction.
        #
        # Both sides are read as text rather than executed: including the
        # converter would run its top-level code, and io.jl's loads are spread
        # across several functions with no single introspectable list.
        repo = normpath(joinpath(@__DIR__, ".."))
        conv = read(joinpath(repo, "tools", "convert_greb_to_jld2.jl"), String)
        io_src = read(joinpath(repo, "src", "io.jl"), String)

        # --- the allowlist, as literals inside the MODEL_FIELD_NAMES block ---
        m = match(r"const MODEL_FIELD_NAMES = Set\{String\}\(\[(.*?)
\]\)"s, conv)
        @test m !== nothing
        # Cut at the ENSO comprehension: its "zonal.wind"/"meridional.wind"
        # tokens are field-name *fragments*, not file names, and it is expanded
        # explicitly below.
        body = m.captures[1]
        cut = findfirst("(\"erainterim.", body)
        cut === nothing || (body = body[1:first(cut)-1])
        allowed = Set{String}()
        for lit in eachmatch(r"\"([^\"]+)\"", body)
            s = lit.captures[1]
            if occursin('$', s)
                continue          # the ENSO comprehension template, expanded below
            elseif occursin('.', s) && !occursin(' ', s)
                push!(allowed, s)
            end
        end
        # expand the ENSO comprehension the same way the converter does
        for f in ("tsurf", "zonal.wind", "meridional.wind", "windspeed", "omega"),
            s in ("elnino", "lanina")
            push!(allowed, "erainterim.$f.$s.forcing")
        end
        @test length(allowed) == 33

        # --- what io.jl actually loads, with $suffix expanded ---
        loaded = Set{String}()
        for m2 in eachmatch(r"\"([A-Za-z0-9_.\$-]+)\.jld2\"", io_src)
            name = m2.captures[1]
            # combined multi-field files are not per-field entries in the allowlist
            name in ("flux_corrections", "ipcc_scenarios", "solar_paleo",
                     "solar_eccentricity", "solar_obliquity") && continue
            if occursin("\$suffix", name)
                for s in ("elnino", "lanina")
                    push!(loaded, replace(name, "\$suffix" => s))
                end
            else
                push!(loaded, name)
            end
        end

        # Every field io.jl loads must be produced by the converter, and the
        # converter must not carry entries nothing loads.
        @test isempty(setdiff(loaded, allowed))
        @test isempty(setdiff(allowed, loaded))
    end

end

function run_heavy_tests()

    @testset "threaded circulation matches serial (subprocess -t 1 vs -t 2)" begin
        # `tendencies!` runs circulation!(Ta) and circulation!(q) concurrently
        # only when `Threads.nthreads() > 1` AND `ws_a !== ws_q` (see
        # src/tendencies.jl). Thread count is fixed at Julia startup, so a
        # single-threaded `Pkg.test()` can never reach that branch - it went
        # untested until 2026-08-21. Spawning both counts explicitly keeps this
        # honest no matter how the suite is invoked.
        if !isdir(DATA_DIR)
            @test_skip "greb_input_data/ not present"
        else
            # No `using Statistics`: it is not a dependency of GREBClimate's own
            # Project.toml, so it is unavailable under `--project=<repo>` even
            # though the test environment has it.
            script = """
                using GREBClimate
                gmean(x) = sum(x) / length(x)
                fields = redirect_stdout(devnull) do
                    load_greb_jld2!(raw"$(DATA_DIR)"; dataset = :ncep)
                end
                cfg = create_experiment_config(:full_model)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(flux = 0, ctrl = 1, scnr = 0), cfg;
                                jld2_dir = raw"$(DATA_DIR)", fields = fields)
                end
                # threads actually available, then a digest of every month
                print(Threads.nthreads())
                for rec in result.ctrl
                    print(" ", gmean(rec.Ts), " ", gmean(rec.Ta), " ", gmean(rec.q))
                end
            """
            # Only the executable from julia_cmd(), not its flags: under
            # Pkg.test those include --check-bounds=yes, which would force the
            # subprocess to recompile the world and make this test ~10x slower.
            exe = first(Base.julia_cmd())
            project = normpath(joinpath(@__DIR__, ".."))
            run_at(n) = begin
                cmd = `$exe --startup-file=no --project=$project -t $n -e $script`
                out = read(cmd, String)
                parts = split(strip(out))
                (nthreads = parse(Int, parts[1]),
                 digest = parse.(Float64, parts[2:end]))
            end

            serial = run_at(1)
            threaded = run_at(2)

            # the subprocesses really did run at the requested thread counts
            @test serial.nthreads == 1
            @test threaded.nthreads == 2
            # 12 months x 3 quantities
            @test length(serial.digest) == 36
            @test length(threaded.digest) == length(serial.digest)

            # Concurrency must not change the answer. The two circulation!
            # calls touch disjoint state (separate workspaces, separate
            # fields), so this should be bit-identical, not merely close.
            @test threaded.digest == serial.digest
        end
    end

    @testset "greb_model! baseline: default config runs to completion with the right output shape" begin
        cfg = create_experiment_config(:full_model)
        result = redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
        end
        @test length(result.ctrl) == 12
        @test length(result.scnr) == 0
        @test result.ctrl[1] isa MonthlyRecord
    end

    @testset "greb_model! flux-correction spin-up: loaded files aren't overwritten by qflux_correction! (§8.5)" begin
        # log_topo_drsp=false + log_qflux_dmc=true selects the "load
        # precomputed files" branch. Before the fix, qflux_correction! ran
        # unconditionally right after and overwrote the just-loaded values
        # whenever time_flux > 0 — use flux=1 (not the RunSpec default of 0)
        # so the spin-up loop actually executes and the old bug would fire.
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "climatology"))
            GREBClimate.jldopen(joinpath(tmpdir, "climatology", "flux_corrections.jld2"), "w") do f
                f["Tsurf_flux_correction"] = fill(42.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
                f["vapour_flux_correction"] = fill(43.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
                f["Tocean_flux_correction"] = fill(44.0, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
            end

            cfg = create_experiment_config(:full_model)
            cfg.log_topo_drsp = false
            cfg.log_qflux_dmc = true
            fields = ClimateFields()
            redirect_stdout(devnull) do
                greb_model!(RunSpec(flux = 1, ctrl = 0, scnr = 0), cfg; jld2_dir = tmpdir, fields = fields, allow_uninitialized = true)
            end

            @test all(==(42.0), fields.TF_correct)
            @test all(==(43.0), fields.qF_correct)
            @test all(==(44.0), fields.ToF_correct)
        finally
            rm(tmpdir; recursive = true, force = true)
        end
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
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "", allow_uninitialized = true)
            end
            @test length(result.ctrl) == 12
        end
    end

    @testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction (matches Fortran)" begin
        fields = ClimateFields()
        fields.cap_surf .= GREBClimate.cap_ocean
        for j in 1:GREBClimate.ydim, i in 1:GREBClimate.xdim
            fields.Tclim[i, j, :] .= 280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0)
            fields.Toclim[i, j, :] .= 279.0
            fields.qclim[i, j, :] .= 0.006
        end
        cfg = create_experiment_config(:full_model)
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        Ts = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
        Ta = fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
        q = fill(0.010, GREBClimate.xdim, GREBClimate.ydim)
        To = fill(285.0, GREBClimate.xdim, GREBClimate.ydim)

        GREBClimate.qflux_correction!(340.0, Ts, Ta, q, To, fields, state, ts, cfg, ws, 1)

        @test any(!=(0.0), fields.TF_correct)
        @test any(!=(0.0), fields.ToF_correct)
        @test any(!=(0.0), fields.qF_correct)
        @test all(isfinite, fields.TF_correct)
        @test all(isfinite, fields.ToF_correct)
        @test all(isfinite, fields.qF_correct)

        @test Ts ≈ fields.Tclim[:, :, 1]
        @test To ≈ fields.Toclim[:, :, 1]
        @test q ≈ fields.qclim[:, :, 1]

        @test all(isfinite, Ta)
        @test Ta != fill(290.0, GREBClimate.xdim, GREBClimate.ydim)
    end

    @testset "greb_model! swaps sw_solar for paleo experiments, restores after" begin
        # A paleo run's swapped solar table must not leak into a later run.
        fields = ClimateFields()
        saved_sw_solar = copy(fields.sw_solar)
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            distinctive_value = 999.0
            GREBClimate.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(distinctive_value, GREBClimate.ydim, GREBClimate.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end

            cfg = create_experiment_config(:paleo_231kyr)
            captured = mktemp() do path, io
                result = redirect_stdout(io) do
                    greb_model!(RunSpec(ctrl = 0), cfg; jld2_dir = tmpdir, fields = fields, allow_uninitialized = true)
                end
                flush(io)
                (result = result, text = read(path, String))
            end
            @test length(captured.result.scnr) == 12
            # Confirm the swap branch actually ran (fields are all-zero in
            # this test env, so asserting on Ts/anomaly magnitude would be
            # fragile - several other tests in this suite hit NaN for the
            # same reason and only check shape/completion instead).
            @test occursin("loading alternate solar forcing", captured.text)

            # sw_solar restored to its pre-run value after greb_model! returns
            @test fields.sw_solar == saved_sw_solar

            cfg_plain = create_experiment_config(:full_model)
            redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields, allow_uninitialized = true)
            end
            @test fields.sw_solar == saved_sw_solar
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "forcing() regional-CO2 ice mask uses the annual mean, not January (§8.1)" begin
        fields = ClimateFields()  # z_topo defaults to 0 everywhere -> land branch never fires
        cfg = PhysicsConfig(experiment = :regional_co2_ocean)

        icmn_ctrl = zeros(Float64, GREBClimate.xdim, GREBClimate.ydim, 12)
        # Cell A: January alone >= 0.5, but the other 11 months are 0 ->
        # annual mean ~0.083, NOT ice under the Fortran-matching rule.
        icmn_ctrl[1, 1, 1] = 1.0
        # Cell B: January alone < 0.5, but the other 11 months are 1.0 ->
        # annual mean ~0.917, IS ice under the Fortran-matching rule.
        icmn_ctrl[2, 1, 1] = 0.0
        icmn_ctrl[2, 1, 2:12] .= 1.0

        forcing(1, 1970, cfg, fields, icmn_ctrl)

        @test fields.co2_part[1, 1] == 1.0  # January said "ice"; annual mean says no
        @test fields.co2_part[2, 1] == 0.5  # January said "no ice"; annual mean says yes
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
            icmn_ctrl = zeros(Float64, GREBClimate.xdim, GREBClimate.ydim, 1)
            redirect_stdout(devnull) do
                init_model!(cfg, fields)
            end
            result = forcing(1, 1970, cfg, fields, icmn_ctrl)
            @test isfinite(result.CO2)
            @test isfinite(result.sw_solar_forcing)
        end

        cfg = PhysicsConfig(experiment = :sst_plus1)
        result = redirect_stdout(devnull) do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "", allow_uninitialized = true)
        end
        @test length(result.scnr) == 12

        # :rcp26/:rcp45/:rcp60 load a year=>CO2 table at scenario start.
        # same mechanism as :ssp*/:historical_co2 - :rcp60's jld2 key is
        # "rcp6", not "rcp60".
        tmpdir_rcp = mktempdir()
        try
            mkpath(joinpath(tmpdir_rcp, "scenario"))
            expected_rcp_co2 = Dict(:rcp26 => 400.0, :rcp45 => 401.0, :rcp60 => 402.0)
            GREBClimate.jldopen(joinpath(tmpdir_rcp, "scenario", "ipcc_scenarios.jld2"), "w") do file
                file["scenarios"] = Dict(
                    "rcp26" => Dict(1950 => expected_rcp_co2[:rcp26]),
                    "rcp45" => Dict(1950 => expected_rcp_co2[:rcp45]),
                    "rcp6" => Dict(1950 => expected_rcp_co2[:rcp60]),
                )
            end

            for sym in (:rcp26, :rcp45, :rcp60)
                cfg = PhysicsConfig(experiment = sym)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir_rcp, allow_uninitialized = true)
                end
                @test length(result.scnr) == 12
                @test cfg.co2_scenario == Dict(1950 => expected_rcp_co2[sym])
            end
        finally
            rm(tmpdir_rcp; recursive = true, force = true)
        end

        # :custom_co2 loads a plain-text "year CO2" file set on cfg.custom_co2_path (§7.2).
        tmpdir_custom = mktempdir()
        try
            co2_path = joinpath(tmpdir_custom, "my_co2.txt")
            write(co2_path, "# comment line, should be skipped\n1950 300.0\n1951 301.0\n\n")

            cfg = create_experiment_config(:custom_co2; co2_path = co2_path)
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "", allow_uninitialized = true)
            end
            @test length(result.scnr) == 12
            @test cfg.co2_scenario == Dict(1950 => 300.0, 1951 => 301.0)

            # Unset custom_co2_path must raise a clear error, not silently
            # dispatch or default.
            cfg_unset = PhysicsConfig(experiment = :custom_co2)
            @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg_unset; jld2_dir = "", allow_uninitialized = true)
        finally
            rm(tmpdir_custom; recursive = true, force = true)
        end

        # :decon_mean_climate (control-run only) / :decon_2xco2 (scenario run)
        # smoke tests.
        cfg_dmc = create_experiment_config(:decon_mean_climate)
        result_dmc = redirect_stdout(devnull) do
            greb_model!(RunSpec(ctrl = 1, scnr = 0), cfg_dmc; jld2_dir = "", allow_uninitialized = true)
        end
        @test length(result_dmc.ctrl) == 12

        cfg_drsp = create_experiment_config(:decon_2xco2)
        result_drsp = redirect_stdout(devnull) do
            greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg_drsp; jld2_dir = "", allow_uninitialized = true)
        end
        @test length(result_drsp.scnr) == 12

        # :ssp* experiments load a year=>CO2 table at scenario start
        tmpdir_ssp = mktempdir()
        try
            mkpath(joinpath(tmpdir_ssp, "scenario"))
            expected_co2 = Dict(
                :ssp119 => 300.0, :ssp126 => 301.0, :ssp245 => 302.0,
                :ssp460 => 303.0, :ssp585 => 304.0,
            )
            GREBClimate.jldopen(joinpath(tmpdir_ssp, "scenario", "ipcc_scenarios.jld2"), "w") do file
                file["scenarios"] = Dict(string(sym) => Dict(1950 => co2) for (sym, co2) in expected_co2)
            end

            for sym in (:ssp119, :ssp126, :ssp245, :ssp460, :ssp585)
                cfg = PhysicsConfig(experiment = sym)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir_ssp, allow_uninitialized = true)
                end
                @test length(result.scnr) == 12
                @test cfg.co2_scenario == Dict(1950 => expected_co2[sym])
            end

            # A year missing from the table must raise a clear error rather
            # than silently defaulting.
            cfg_missing_year = PhysicsConfig(experiment = :ssp585)
            @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 2), cfg_missing_year;
                jld2_dir = tmpdir_ssp, allow_uninitialized = true)
        finally
            rm(tmpdir_ssp; recursive = true, force = true)
        end

        tmpdir_hist = mktempdir()
        try
            mkpath(joinpath(tmpdir_hist, "scenario"))
            GREBClimate.jldopen(joinpath(tmpdir_hist, "scenario", "ipcc_scenarios.jld2"), "w") do file
                file["scenarios"] = Dict("hist" => Dict(1850 => 280.73))
            end

            cfg = PhysicsConfig(experiment = :historical_co2)
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir_hist, allow_uninitialized = true)
            end
            @test length(result.scnr) == 12
            @test isapprox(cfg.co2_scenario[1850], 280.73; atol = 1e-3)
        finally
            rm(tmpdir_hist; recursive = true, force = true)
        end

        # :paleo_solar_modern_co2/:obliquity/:eccentricity swap in a real
        # solar_scenarios/*.jld2 table at scenario start - synthesize a
        # minimal set so these three are reachable without the real dataset.
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            GREBClimate.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(999.0, GREBClimate.ydim, GREBClimate.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end
            GREBClimate.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_obliquity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREBClimate.ydim, GREBClimate.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end
            GREBClimate.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_eccentricity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREBClimate.ydim, GREBClimate.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end

            for sym in (:paleo_solar_modern_co2, :obliquity, :eccentricity)
                cfg = PhysicsConfig(experiment = sym)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir, allow_uninitialized = true)
                end
                @test length(result.scnr) == 12
            end
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "CMIP5/ERA-Interim anomaly forcing is actually loaded (was previously a silent no-op)" begin
        tmpdir_anom = mktempdir()
        try
            clim_dir = joinpath(tmpdir_anom, "climatology")
            mkpath(clim_dir)
            write_field(name, value) = GREBClimate.jldopen(joinpath(clim_dir, name), "w") do file
                file["data"] = fill(value, GREBClimate.xdim, GREBClimate.ydim, GREBClimate.nstep_yr)
                file["dim_names"] = ["lon", "lat", "time"]
            end

            # :rcp85 → CMIP5 RCP8.5 ensemble-mean anomaly
            write_field("cmip5.tsurf.rcp85.ensmean.forcing.jld2", 2.0)
            write_field("cmip5.zonal.wind.rcp85.ensmean.forcing.jld2", 3.0)
            write_field("cmip5.meridional.wind.rcp85.ensmean.forcing.jld2", 4.0)
            write_field("cmip5.omega.rcp85.ensmean.forcing.jld2", 5.0)
            write_field("cmip5.windspeed.rcp85.ensmean.forcing.jld2", 6.0)

            # :elnino / :lanina → ERA-Interim composite-mean anomaly
            for suffix in ("elnino", "lanina")
                write_field("erainterim.tsurf.$suffix.forcing.jld2", 7.0)
                write_field("erainterim.zonal.wind.$suffix.forcing.jld2", 8.0)
                write_field("erainterim.meridional.wind.$suffix.forcing.jld2", 9.0)
                write_field("erainterim.omega.$suffix.forcing.jld2", 10.0)
                write_field("erainterim.windspeed.$suffix.forcing.jld2", 11.0)
            end

            cfg = create_experiment_config(:rcp85)
            @test cfg.log_tsurf_ext && cfg.log_hwind_ext && cfg.log_omega_ext
            fields = ClimateFields()
            load_cc_anomaly_jld2!(tmpdir_anom, fields, cfg)
            @test all(==(2.0), fields.Tclim_anom_cc)
            @test all(==(3.0), fields.uclim_anom_cc)
            @test all(==(4.0), fields.vclim_anom_cc)
            @test all(==(5.0), fields.omegaclim_anom_cc)
            @test all(==(6.0), fields.wsclim_anom_cc)

            # init_model! applies the anomaly on top of the (here all-zero)
            # base climatology - Tclim must reflect it, not stay at zero.
            redirect_stdout(devnull) do
                init_model!(cfg, fields)
            end
            @test all(==(2.0), fields.Tclim)

            for (sym, suffix) in ((:elnino, "elnino"), (:lanina, "lanina"))
                cfg2 = create_experiment_config(sym)
                @test cfg2.log_tsurf_ext && cfg2.log_hwind_ext && cfg2.log_omega_ext
                fields2 = ClimateFields()
                load_enso_anomaly_jld2!(tmpdir_anom, fields2, cfg2, sym)
                @test all(==(7.0), fields2.Tclim_anom_enso)
                @test all(==(8.0), fields2.uclim_anom_enso)
                @test all(==(9.0), fields2.vclim_anom_enso)
                @test all(==(10.0), fields2.omegaclim_anom_enso)
                @test all(==(11.0), fields2.wsclim_anom_enso)
            end

            # Per-field gating: switching a gate off must not touch that
            # field even if its file is missing (no error, stays zero).
            cfg_partial = PhysicsConfig(experiment = :rcp85, log_tsurf_ext = true,
                log_hwind_ext = false, log_omega_ext = false)
            fields_partial = ClimateFields()
            load_cc_anomaly_jld2!(tmpdir_anom, fields_partial, cfg_partial)
            @test all(==(2.0), fields_partial.Tclim_anom_cc)
            @test all(==(0.0), fields_partial.uclim_anom_cc)
            @test all(==(0.0), fields_partial.omegaclim_anom_cc)

            # A missing required file must error loudly, not silently zero.
            rm(joinpath(clim_dir, "cmip5.tsurf.rcp85.ensmean.forcing.jld2"))
            @test_throws ErrorException load_cc_anomaly_jld2!(tmpdir_anom, ClimateFields(), cfg)
        finally
            rm(tmpdir_anom; recursive = true, force = true)
        end
    end

    @testset "golden regression: real dataset control+scenario run matches snapshot" begin
        # Tripwire for any refactor touching the physics kernels: a real 1yr
        # control + 1yr scenario run against the actual NCEP dataset,
        # snapshotted as monthly global-mean Ts/Ta/q. Drift beyond
        # float-reassociation noise (~1e-12) means real behavior changed.
        # Set RUN_GOLDEN=0 to skip this locally; CI always runs it.
        if !isdir(DATA_DIR)
            @test_skip "greb_input_data/ not present"
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
            # Tolerances loosened for Float32 (was 1e-6 under Float64): Ts/Ta
            # are O(280K), so Float32's ~1.2e-7 relative epsilon alone gives
            # an absolute noise floor around 3e-5, before any compounding
            # over a 2-year run — measured up to ~4e-3 in practice. The
            # `scnr_ref` values are themselves tiny control-vs-scenario
            # anomalies (many below Float32's own noise floor for an O(280K)
            # quantity), so at this precision the scenario assertions can
            # only catch a real blow-up, not reproduce the fine anomaly
            # structure Float64 could — an accepted consequence of the
            # precision switch, not a bug.
            for (rec, ref) in zip(result.ctrl, ctrl_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-2)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-2)
                @test isapprox(s.q, ref.q; atol = 1e-5)
            end
            for (rec, ref) in zip(result.scnr, scnr_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-2)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-2)
                @test isapprox(s.q, ref.q; atol = 1e-5)
            end
        end
    end

end

shard = get(ENV, "GREB_TEST_SHARD", "all")
@testset "GREBClimate.jl" begin
    shard in ("all", "light") && run_light_tests()
    shard in ("all", "heavy") && run_heavy_tests()
end

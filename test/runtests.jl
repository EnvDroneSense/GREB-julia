using GREB
using Test

# Smoke tests that do NOT require the (large, external) JLD2 input data.
# They check that the package loads, its types build, and the grid/constants
# are intact after the notebook -> package extraction. Full integration runs
# (which need `greb_dataset_jld2/`) are demonstrated in examples/run_greb.jl.

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

    @testset "set_hydrology_parameters! writes into cfg, not a disconnected global" begin
        # Regression: this used to declare `global c_q, c_rq, c_omega, c_omegastd`
        # and assign those instead of cfg.c_q etc. hydro! reads cfg.c_q, so
        # log_rain silently had zero effect on every run regardless of its
        # value. Check cfg's own fields match the HYDRO_PARAMS table directly.
        for log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_rain = log_rain
            set_hydrology_parameters!(cfg)
            expected = GREB.HYDRO_PARAMS[log_rain]
            @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == expected
        end

        # NCEP adjustment branch (log_rain == 0 && log_clim == 1) overrides the table
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 0
        cfg.log_clim = 1
        set_hydrology_parameters!(cfg)
        @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == (-1.27, 1.99, -16.54, 21.15)
    end

    @testset "workspaces & accumulators" begin
        ws = CirculationWorkspace()
        @test ws isa CirculationWorkspace

        acc = MonthlyAccumulator()
        @test acc isa MonthlyAccumulator
        @test (GREB.reset!(acc); true)   # reset! runs without error

        ts = TimeState(1, 1)
        @test ts.jday == 1
        @test ts.ityr == 1
    end

    @testset "MonthlyRecord type" begin
        @test MonthlyRecord <: NamedTuple
        @test :Ts in fieldnames(MonthlyRecord)
        @test :precip in fieldnames(MonthlyRecord)
    end

    @testset "greb_model! runs without notebook globals" begin
        # Regression: qflux_correction!/greb_model! used to reference the Pluto
        # @bind globals `time_flux`/`jld2_dir`. They are now parameters, so the
        # model must run to completion on default (unloaded) fields. Values are
        # NaN without real JLD2 data — we only assert it runs and shapes are OK.
        cfg = create_experiment_config(:full_model)
        result = redirect_stdout(devnull) do
            greb_model!(0, 1, 0, cfg; jld2_dir = "")
        end
        @test length(result.ctrl) == 12
        @test length(result.scnr) == 0
        @test result.ctrl[1] isa MonthlyRecord
    end

    @testset "greb_model! runs across log_eva / log_rain branches" begin
        # Regression for a `ws.cE` vs `ws.cE_buf` field-name typo in hydro!'s
        # log_eva == 0 branch: it threw `type CirculationWorkspace has no
        # field cE` at runtime, but was invisible to the test above because
        # PhysicsConfig's default (log_eva = -1) never reached that branch.
        for log_eva in (-1, 0, 1), log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            cfg.log_rain = log_rain
            result = redirect_stdout(devnull) do
                greb_model!(0, 1, 0, cfg; jld2_dir = "")
            end
            @test length(result.ctrl) == 12
        end
    end

    @testset "read_jld2 rejects non-JLD2 input" begin
        tmp = tempname() * ".jld2"
        write(tmp, "not a jld2 file")
        @test_throws Exception read_jld2(tmp)
        rm(tmp; force = true)
    end

end

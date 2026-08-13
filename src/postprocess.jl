"""
    build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}

Returns a 12-month climatology taken from the *final* year of `records`.
"""
function build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(records) && return MonthlyRecord[]

    nmonths = 12
    n = length(records)
    final_year = records[max(1, n - nmonths + 1):n]

    clim = MonthlyRecord[]
    for mon in 1:nmonths
        push!(clim, mon <= length(final_year) ? final_year[mon] : final_year[1])
    end
    return clim
end

"""
    apply_scenario_anomalies(scnr_records, ctrl_clim)::Vector{MonthlyRecord}

Subtracts the matching calendar month of `ctrl_clim` (from
[`build_monthly_climatology`](@ref)) from each record in `scnr_records`,
turning absolute monthly output into anomalies relative to the control run.
"""
function apply_scenario_anomalies(scnr_records::Vector{MonthlyRecord}, ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(scnr_records) && return scnr_records
    isempty(ctrl_clim) && return scnr_records

    fields = propertynames(scnr_records[1])
    anom = MonthlyRecord[]

    for (idx, rec) in enumerate(scnr_records)
        mon = mod(idx - 1, 12) + 1
        ref = ctrl_clim[mon]
        push!(anom, NamedTuple{fields}(
            Tuple(getfield(rec, fld) .- getfield(ref, fld) for fld in fields)))
    end
    return anom
end

"""
    compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})

Returns `ctrl_output`'s `ice` field from the *final* year, as an
`(xdim, ydim, 12)` array.
"""
function compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float32, xdim, ydim, 12)
    isempty(ctrl_output) && return ice_months

    nmonths = 12
    n = length(ctrl_output)
    final_year = ctrl_output[max(1, n - nmonths + 1):n]

    for (mon, rec) in enumerate(final_year)
        @. ice_months[:, :, mon] = rec.ice   # note: field name is `ice`, not `ice_cover`
    end

    return ice_months
end

const df = let
    Fc = 12
    n = 2
    bw = Butterworth(n)
    corrfac = inv((2^inv(n)-1)^(1/4)) # Correction factor for Fc of multi-pass filters
    lpf = Lowpass(Fc * corrfac; fs=100)
    digitalfilter(lpf, bw)
end

const lab_orientation = [ 0  0 -1;
                         -1  0  0;
                          0  1  0 ]

function getheel(trial::Trial; kwargs...)
    getheel(readsource(getsource(trial, Source{C3DFile})); kwargs...)
end

function getheel(file::C3DFile;
    lheelmkr="LHEE",
    rheelmkr="RHEE",
    axis=Colon(),
    start=firstindex(file.point[lheelmkr],1),
    finish=lastindex(file.point[lheelmkr],1)
)
    lheel = file.point[lheelmkr]
    rheel = file.point[rheelmkr]

    frheel = filtfilt(df, rheel[start:finish,axis])
    flheel = filtfilt(df, lheel[start:finish,axis])

    return flheel, frheel
end

function roerdink2008(trial; fo_minprom=1200)
    c3dsrc = readsource(getsource(trial, Source{C3DFile}); strip_prefixes=true)
    fs = c3dsrc.groups[:POINT][Int, :RATE]

    lheel, rheel = getheel(c3dsrc; axis=3)
    lfcpred, _ = peakproms!(argminima(lheel, 10), lheel; minprom=30)
    rfcpred, _ = peakproms!(argminima(rheel, 10), rheel; minprom=30)

    lheel_vel = centraldiff(lheel; dt=inv(fs), padding=ForwardBackwardPad())
    rheel_vel = centraldiff(rheel; dt=inv(fs), padding=ForwardBackwardPad())

    lfopred, _ = peakproms!(argmaxima(lheel_vel, 10), lheel_vel; minprom=fo_minprom)
    rfopred, _ = peakproms!(argmaxima(rheel_vel, 10), rheel_vel; minprom=fo_minprom)

    return Dict("LFC" => totimes(lfcpred, fs), "RFC" => totimes(rfcpred, fs),
                "LFO" => totimes(lfopred, fs), "RFO" => totimes(rfopred, fs))
end

notnanmissing(x) = isnan(x) === false

function validsegments(x::AbstractVector)
    segs = UnitRange{Int}[]
    bi = li = firstindex(x)
    lv = x[li]

    for i in eachindex(x)
        if notnanmissing(lv)
            if !notnanmissing(x[i])
                push!(segs, bi:li)
            end
        elseif notnanmissing(x[i])
            bi = i
        end
        li = i
        lv = x[li]
    end
    if notnanmissing(lv)
        push!(segs, bi:li)
    end

    return segs
end

filterpoint_zerolag(p::AbstractArray, f) = _filterpoint_zerolag!(copy(p), p, f)
filterpoint_zerolag!(p::AbstractArray, f) = _filterpoint_zerolag!(p, copy(p), f)

function _filterpoint_zerolag!(pf, p, f; minlength=50)
    segs = validsegments(@view(p[:,1]))
    filter!(x -> length(x) > minlength, segs)

    for seg in segs
        @views pf[seg, :] = filtfilt(f, p[seg, :])
    end

    return pf
end

function filterpoints!(file::C3DFile, f)
    filterpoint_zerolag!.(values(file.point), Ref(f))

    return file
end


using HypothesisTests, Printf, Distributions

struct TwoOneSidedTTest <: HypothesisTests.TTest
    main::OneSampleTTest
    upper::OneSampleTTest
    lower::OneSampleTTest
    pooledvar::Float64

end

Base.getproperty(x::TwoOneSidedTTest, p::Symbol) = p ∈ fieldnames(OneSampleTTest) ? getfield(x.main, p) : getfield(x, p)


function TwoOneSidedTTest(x, y, Δ=0; d=nothing, ΔL=-Δ, ΔU=Δ)
    if !isnothing(d)
        diff_sd = std(x .- y)
        ΔL = -d*diff_sd
        ΔU = d*diff_sd
    end
    main = OneSampleTTest(x,y)
    upper = OneSampleTTest(x,y,ΔU)
    lower = OneSampleTTest(x,y,ΔL)

    pvar = (var(x) + var(y))/2

    return TwoOneSidedTTest(main, upper, lower, pvar)
end

HypothesisTests.testname(::TwoOneSidedTTest) = "TOST Paired Samples t-test"
HypothesisTests.population_param_of_interest(x::TwoOneSidedTTest) = HypothesisTests.population_param_of_interest(x.main)

HypothesisTests.confint(x::TwoOneSidedTTest; level=0.95, tail=:both) = confint(x.main; level, tail)
HypothesisTests.pvalue(x::TwoOneSidedTTest; tail=:both) = pvalue(x.main; tail)

function tost_pvalue(x::TwoOneSidedTTest)
    return max(pvalue(x.lower; tail=:right), pvalue(x.upper; tail=:left))
end

function LEAD(x::TwoOneSidedTTest; α=0.05)
    ad = cquantile(TDist(x.main.n*2-2), α)*√(x.pooledvar*2/x.main.n)
    lead_lower = x.main.xbar + ad
    lead_upper = -x.main.xbar + ad

    return max(abs(lead_lower), abs(lead_upper))
end

function HypothesisTests.show_params(io::IO, x::TwoOneSidedTTest, indent="")
    println(io, indent, "number of observations:   $(x.main.n)")
    println(io, indent, "t-statistic:              $(x.main.t)")
    println(io, indent, "degrees of freedom:       $(x.main.df)")
    println(io, indent, "empirical standard error: $(x.main.stderr)")
    println(io, indent, "equivalence bounds:       ΔL: $(x.lower.μ0), ΔU: $(x.upper.μ0)")

    lowerp = pvalue(x.lower; tail=:right)
    upperp = pvalue(x.upper; tail=:left)
    @printf io "\n%sTOST outcome with 95%% confidence: %sreject TOST h_0\n" indent ((lowerp ≤ .05 && upperp ≤ .05) ? "" : "fail to ")
    @printf io "%sTOST Lower:               t(%d)=%.2f, p%s%s\n" indent x.lower.df x.lower.t (lowerp < .001 ? '<' : '=') lstrip(string(round(clamp(lowerp, .001, 1.); digits=3)), '0')
    @printf io "%sTOST Upper:               t(%d)=%.2f, p%s%s\n" indent x.upper.df x.upper.t (upperp < .001 ? '<' : '=') lstrip(string(round(clamp(upperp, .001, 1.); digits=3)), '0')
    # println(io, indent, "LEAD: [max(|cₗ|, |cᵤ|)]    $(LEAD(x))")

    return nothing
end

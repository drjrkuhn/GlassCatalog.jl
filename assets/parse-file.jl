using DataArrays, DataFrames, YAML, Interpolations

@enum(INDEXFN,
  Sellmeier =1,
  Sellmeier2 =2,
  Polynomial =3 ,
  RFInfo =4,
  Cauchy =5,
  Gases =6,
  Herzberger =7,
  Retro =8,
  Exotic =9
)

const λ_i = (0.36501)	#	Hg	UV
const λ_h = (0.40466)	#	Hg	violet
const λ_g = (0.43584)	#	Hg	blue
const λ_Fp = (0.47999)	#	Cd	blue
const λ_F = (0.48613)	#	H	blue
const λ_e = (0.54607)	#	Hg	green
const λ_d = (0.58756)	#	He	yellow
const λ_D = (0.5893)	#	Na	yellow
const λ_Cp = (0.64385)	#	Cd	red
const λ_C = (0.65627)	#	H	red
const λ_r = (0.70652)	#	He	red
const λ_Ap = (0.7682)	#	K	IR
const λ_s = (0.85211)	#	Cs	IR
const λ_t = (0.101398)	#	Hg	IR

_fraunhofer_lines = Dict(
"i" => ("Hg", 0.36501, "blue"),
"h" => ("Hg", 0.40466, "blue"),
"g" => ("Hg", 0.43585, "blue"),
"F'" => ("Cd", 0.47999, "green"),
"F" => ("Hβ", 0.486134, "infrared"),
"e" => ("Hg", 0.546073, "infrared"),
"d" => ("He", 0.58756, "infrared"),
"D" => ("Na", 0.5893, "red"),
"C'" => ("Cd", 0.64385, "red"),
"C" => ("Hα", 0.656281, "red"),
"r" => ("He", 0.70652, "ultraviolet"),
"A'" => ("K", 0.7682, "violet"),
"s" => ("Cs", 0.85211, "yellow"),
"t" => ("Hg", 1.01398, "yellow"),
"t2" => ("Ni", 0.299444, ""),
"T" => ("Fe", 0.302108, ""),
"P" => ("Ti+", 0.336112, ""),
"N" => ("Fe", 0.358121, ""),
"L" => ("Fe", 0.382044, ""),
"K" => ("Ca+", 0.3933666, ""),
"H" => ("Ca+", 0.396847, ""),
"h2" => ("Hδ", 0.410175, ""),
"G2" => ("Ca", 0.430774, ""),
"G" => ("Fe", 0.43079, ""),
"G'" => ("Hγ", 0.434047, ""),
"e2" => ("Fe", 0.438355, ""),
"d2" => ("Fe", 0.466814, ""),
"c" => ("Fe", 0.495761, ""),
"b4" => ("Mg", 0.516733, ""),
"b3" => ("Fe", 0.516891, ""),
"b2" => ("Mg", 0.51727, ""),
"b1" => ("Mg", 0.518362, ""),
"E2" => ("Fe", 0.527039, ""),
"D3" => ("He", 0.5875618, ""),
"D2" => ("Na", 0.588995, ""),
"D1" => ("Na", 0.589592, ""),
"a" => ("O2", 0.627661, ""),
"B" => ("O2", 0.686719, ""),
"A" => ("O2", 0.75937, ""),
"Z" => ("O2", 0.822696, ""),
"y" => ("O2", 0.898765, "")
)

fraunhofer_λum(line::AbstractString) = get(_fraunhofer_lines,line,NaN)[2]
fraunhofer_source(line::AbstractString) = get(_fraunhofer_lines,line,"")[1]
fraunhofer_color(line::AbstractString) = get(_fraunhofer_lines,line,"")[3]
fraunhofer_standards() = keys(filter((k,v)->!isempty(v[3]), _fraunhofer_lines))

function refindex{T<:Real}(::Type{Val{Sellmeier}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local λsq::T = λum^2
  local nsq::T = one(T) + C[1]
  for i = 2:2:N-1
    nsq += C[i] * λsq / (λsq - C[i+1]^2)
  end
  sqrt(nsq)
end

function refindex{T<:Real}(::Type{Val{Sellmeier2}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local λsq::T = λum^2
  local nsq::T = one(T) + C[1]
  for i = 2:2:N-1
    nsq += C[i] * λsq / (λsq - C[i+1])
  end
  sqrt(nsq)
end

function refindex{T<:Real}(::Type{Val{Polynomial}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local nsq::T = C[1]
  for i = 2:2:N-1
    nsq += C[i] * λum^C[i+1]
  end
  sqrt(nsq)
end

function refindex{T<:Real}(::Type{Val{RFInfo}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=11 "index equation requires an odd number of coefficients ≥ 11"
  local λsq::T = λum^2
  local nsq::T = C[1]
  for i = 2:4:9
    nsq += C[i] * λum^C[i+2] / (λsq - C[i+3]^C[i+4])
  end
  for i = 10:2:N-1
    nsq += C[i]*λum^C[i+1]
  end
  sqrt(nsq)
end

function refindex{T<:Real}(::Type{Val{Cauchy}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local n::T = C[1]
  for i = 2:2:N-1
    n += C[i]*λum^C[i+1]
  end
  n
end

function refindex{T<:Real}(::Type{Val{Gases}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local n::T = one(T) + C[1]
  local invsqrtλ::T = λum^-2
  for i = 2:2:N-1
    n += C[i] / (C[i+1] - invsqrtλ)
  end
  n
end

function refindex{T<:Real}(::Type{Val{Herzberger}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert N>=6 "index equation requires 6 coefficients"
  local n::T = C[1]
  local λsq::T = λum^2
  n += C[2] / (λsq - 0.028)
  n += C[3]*(1/(λsq - 0.028))^2
  for i = 4::N
    n += C[i] * λum^(2*(i-3))
  end
  n
end

function refindex{T<:Real}(::Type{Val{Retro}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert N>=4 "index equation requires 4 coefficients"
  local λsq::T = λum^2
  local r::T = C[1] + C[2]*λsq/(λsq-C[3]) + C[4]*λsq
  sqrt(-2r - 1) / sqrt(r - 1)
end

function refindex{T<:Real}(::Type{Val{Exotic}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert N>=6 "index equation requires 6 coefficients"
  local λsq::T = λum^2
  local nsq::T = C[1] + C[2]/(λsq-C[3]) + C[4]*(λum - C[5])/((λum - C[5])^2 + C[6])
  sqrt(nsq)
end


fusedSiCstr = "0 0.6961663 0.0684043 0.4079426 0.1162414 0.8974794 9.896161"
fusedSiC = readdlm(IOBuffer(fusedSiCstr),' ', Float64)

type Material
  name::AbstractString
  nD::Float64
  VD::Float64
  fnum::Int
  range::Tuple{Float64, Float64}
  C::Array{Float64,1}
  Ntable::Array{Float64,2}
  Ktable::Array{Float64,2}
  Ninterp::Nullable{Interpolations.GriddedInterpolation}
  Kinterp::Nullable{Interpolations.GriddedInterpolation}
  info::Dict{AbstractString,Any}
  comments::AbstractString
  references::AbstractString
  Material(name) = new(name)
end

function Material(name::AbstractString, path::AbstractString, file::AbstractString)
  matd = YAML.load_file(joinpath(path,file))
  m = Material(name)
  m.name = name
  m.nD = NaN
  m.VD = NaN
  m.fnum = -1
  m.range = (NaN, NaN,)
  m.C = Array(Float64,0)
  m.Ntable = Array(Float64,0,0)
  m.Ktable = Array(Float64,0,0)
  m.Ninterp = Nullable()
  m.Kinterp = Nullable()
  if haskey(matd, "DATA")
    for data in matd["DATA"]
      typestr = lowercase(get(data, "type", ""))
      if startswith(typestr,"formula") && haskey(data, "coefficients")
        m.fnum = parse(Int,split(typestr)[end])
        m.C = vec(readdlm(IOBuffer(data["coefficients"]), ' ', Float64)[1,:])
        range = haskey(data, "range") ? readdlm(IOBuffer(data["range"]),' ', Float64) : [NaN, NaN]
        m.range = tuple(range...)
      elseif startswith(typestr,"tabulated") && haskey(data,"data")
        tabnk = lowercase(split(typestr)[end])
        hasn = search(tabnk, 'n') > 0
        hask = search(tabnk, 'k') > 0
        tab = readdlm(IOBuffer(data["data"]), ' ', Float64)
        if any(isnan, m.range)
          m.range = extrema(tab[:,1])
        end
        if hasn && hask
          m.Ntable = tab[:,[1,2]]
          m.Ktable = tab[:,[1,3]]
        elseif hasn
          m.Ntable = tab
        else
          m.Ktable = tab
        end
      end
    end
  end
  if length(m.Ntable) > 0
    knots = (vec(m.Ntable[:,1]),)  # tuple of x coords
    vals = vec(m.Ntable[:,2])      # f(x)
    m.Ninterp = Nullable(interpolate(knots, vals, Gridded(Linear())))
  end
  if length(m.Ktable) > 0
    knots = (vec(m.Ktable[:,1]),)  # tuple of x coords
    vals = vec(m.Ktable[:,2])      # f(x)
    m.Kinterp = Nullable(interpolate(knots, vals, Gridded(Linear())))
  end

  m.info = get(matd, "INFO", Dict())
  m.comments = get(matd, "COMMENTS", "")
  m.references = get(matd, "REFERENCES","")
  m.nD = get(m.info, "n_d", NaN)
  m.VD = get(m.info, "V_d", NaN)
  #println(">> range: $(m.range[1]) to $(m.range[2])")
  if (isnan(m.nD) || isnan(m.VD)) && m.range[1] <= λ_F && λ_C <= m.range[2]
    nD = refindex(m, λ_D)
    nC = refindex(m, λ_C)
    nF = refindex(m, λ_F)
    #nF = refindex(Val{INDEXFN(m.fnum)}, m.C, λ_F)
    #println(">>Calculated nD=$nD, nC=$nC, nF=$nF")
    m.nD = nD
    m.VD = (m.nD - 1) / (nF - nC)
  end
  m
end

function refindex{T<:Real}(m::Material, λum::T)
  if 1 <= m.fnum <= 9 && length(m.C) > 0
    refindex(Val{INDEXFN(m.fnum)}, m.C, λum)
  elseif !isnull(m.Ninterp) && m.range[1] <= λum <= m.range[2]
    get(m.Ninterp)[λum]
  else
    NaN
  end
end


dbfolder = Pkg.dir("GlassCatalog","assets","refractiveindex.info-database","database")

println("=== L-LAH86 ===")
lah86 = Material("L-LAH86", dbfolder, "glass/ohara/L-LAH86.yml")
dump(lah86)
println()

println("=== BSC7 ===")
bsc7 = Material("BSC7", dbfolder, "glass/hoya/BSC7.yml")
dump(bsc7)
println()

println("=== GaInP ===")
gainp = Material("GaInP", dbfolder, "other/mixed crystals/GaInP/Schubert.yml")
dump(gainp)
println()

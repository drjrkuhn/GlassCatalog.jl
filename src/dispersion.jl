module dispersion

using Interpolations

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

line_λum(line::AbstractString) = get(_fraunhofer_lines,line,NaN)[2]
line_source(line::AbstractString) = get(_fraunhofer_lines,line,"")[1]
line_color(line::AbstractString) = get(_fraunhofer_lines,line,"")[3]
line_standards() = keys(filter((k,v)->!isempty(v[3]), _fraunhofer_lines))

# standard dispersion calculation wavelengths
const λum_D = line_λum("D")
const λum_C = line_λum("C")
const λum_F = line_λum("F")

# forumulas used by the database
@enum(INDEXFN, Sellmeier =1, Sellmeier2 =2, Polynomial =3, RFInfo =4,
  Cauchy =5, Gases =6, Herzberger =7, Retro =8, Exotic =9)

# Sellmeier (formula 1)
function refindex{T<:Real}(::Type{Val{Sellmeier}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local λsq::T = λum^2
  local nsq::T = one(T) + C[1]
  for i = 2:2:N-1
    nsq += C[i] * λsq / (λsq - C[i+1]^2)
  end
  return sqrt(nsq)
end

# Sellmeier2 (formula 2)
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

# Polynomial (formula 3)
function refindex{T<:Real}(::Type{Val{Polynomial}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local nsq::T = C[1]
  for i = 2:2:N-1
    nsq += C[i] * λum^C[i+1]
  end
  sqrt(nsq)
end

# RFInfo (formula 4)
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

# Cauchy (formula 5)
function refindex{T<:Real}(::Type{Val{Cauchy}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert isodd(N) && N>=3 "index equation requires an odd number of coefficients"
  local n::T = C[1]
  for i = 2:2:N-1
    n += C[i]*λum^C[i+1]
  end
  n
end

# Gasses (formula 6)
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

# Herzberger (formula 7)
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

# Retro (formula 8)
function refindex{T<:Real}(::Type{Val{Retro}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert N>=4 "index equation requires 4 coefficients"
  local λsq::T = λum^2
  local r::T = C[1] + C[2]*λsq/(λsq-C[3]) + C[4]*λsq
  sqrt(-2r - 1) / sqrt(r - 1)
end

# Exotic (formula 9)
function refindex{T<:Real}(::Type{Val{Exotic}}, C::Array{T,1}, λum::T)
  local N = length(C)
  @assert N>=6 "index equation requires 6 coefficients"
  local λsq::T = λum^2
  local nsq::T = C[1] + C[2]/(λsq-C[3]) + C[4]*(λum - C[5])/((λum - C[5])^2 + C[6])
  sqrt(nsq)
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

immutable SeriesInterp{T<:Real}
  extrema::Tuple{T, T}
  interp::Interpolations.GriddedInterpolation
  function SeriesInterp(xvals::Array{T,1}, yvals::Array{T,1})
    new(extrema(xvals), interpolate((xvals,), vals, Gridded(Linear())))
  end
  SeriesInterp(table::Array{T, 2}) = SeriesInterp(vec(table[:,1]), vec(table[:,2]))
end

interpolate{T<:Real}(intp::SeriesInterp, x::T) = extrema[1] <= x <= extrema[2] ? intp.interp[x] : NaN


type Material{T<:Real}
  name::AbstractString
  nD::T
  VD::T
  fnum::Int
  extrema::Tuple{T, T}
  C::Array{T,1}
  nintp::Nullable{SeriesInterp{T}}
  kintp::Nullable{SeriesInterp{T}}
  info::Dict{AbstractString,Any}
  comments::AbstractString
  references::AbstractString
  Material(name) = new(name)
end


end # module

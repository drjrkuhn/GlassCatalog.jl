
abstract IndexFormula{T<:Number}

"""
refindex(formula, λ_μm) computes the refractive index for a given wavelength in μm.
"""
refindex{T}(::IndexFormula{T}, λ_μm::T) = zero(T);


type SellmeierFormula{T<:Number} <: IndexFormula{T}
  n::Int8
  A::T
  B::Array{T,1}
  C::Array{T,1}
  λn::Bool
  function SellmeierFormula(A,B,C,λn=true)
    n = length(B)
    if length(C) != n
      error("Coefficients B and C of different length")
    else
      new(n,A,B,C,λn)
    end
  end
end


function refindex{T}(formula::SellmeierFormula{T}, λ_μm::T)
  nsq=one(T)+formula.A
  λsq = λ_μm^2
  for i = 1:formula.n
    num = i==formula.n && !formula.λn ? formula.B[i] : formula.B[i] * λsq
    nsq += num/(λsq - formula.C[i])
  end
  return sqrt(nsq)
end

bk7 = SellmeierFormula{Float64}(0, [1.03961212, 0.231792344, 1.01046945], [6.00069867e-3, 2.00179144e-2, 1.03560653e2])
refindex(bk7,0.5)

using Gadfly
plot(x->refindex(bk7,x),0.21,1.6)

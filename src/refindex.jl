type DispersionFormula
  ncoef::Integer
end

type SellmeierFormula <: DispersionFormula
  SellmeierFormula() = super(17)
end

type Sellmeier2Formula <: DispersionFormula
  Sellmeier2Formula() = super(17)
end

type PolynomialFormula <: DispersionFormula
  PolynomialFormula() = super(17)
end

type CauchyFormula <: DispersionFormula
  CauchyFormula() = super(11)
end

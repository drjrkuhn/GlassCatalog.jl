module Test

using DataFrames

df = DataFrame(a=ASCIIString[], b=Int[])

push!(df, ["Hello", 1])
push!(df, ["World", 2])

show(df)

end

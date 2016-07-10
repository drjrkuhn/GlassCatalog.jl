using Taro, DataFrames

file = joinpath(Pkg.dir("GlassCatalog"),"assets","schott","schott_optical_glass_catalogue_excel_may_2016.xls");

firstrow = 4
lastrow = 20
firstcol = "A"
lastcol = "FJ"
range = "$firstcol$firstrow:$lastcol$lastrow";

table = Taro.readxl(file,"Datatable",range)

typeof(table)

#describe(table)

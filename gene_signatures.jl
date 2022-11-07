# module Init
using Pkg
if isfile("Project.toml") 
    if !isfile("Manifest.toml")
        Pkg.instantiate()
    else 
        Pkg.activate(".")
    end 
end
using CSV
using BSON
using DataFrames

data = CSV.read("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", DataFrame, header = 3, delim = "\t")
CSV.read(dsa)


struct Data 
    data::Array{Float64, 2}
    cols::Array{String31, 1}
    rows::Array{String31, 1}
end 

values = Matrix(data[:,3:end])
GTEX_tpm = Data(values, names(data)[3:end], data[:, "Description"])

GTEX_tpm_test = Data(values[:,10000:end], names(data)[10003:end], data[:, "Description"])
bson("data/gtex_tpm.bson", Dict("data" => GTEX_tpm_test))
for i in 1:size(values)[1]
    if typeof(values[1,:]) !=  Vector{Float64}
        println(i)
    end 
end
typeof(values)
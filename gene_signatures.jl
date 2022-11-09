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
using Statistics
using HDF5
struct Data 
    data::Array{Float64, 2}
    cols::Array{String31, 1}
    rows::Array{String31, 1}
end 
data = CSV.read("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", DataFrame, header = 3, delim = "\t")
gencode = CSV.read("/u/leucegene/data/Homo_sapiens.GRCh38_H32/annotations/Homo_sapiens.GRCh38.Gencode37.genes.tsv",DataFrame,  delim = "\t")
ensmbl = CSV.read("/u/leucegene/data/Homo_sapiens.GRCh38_H32/annotations/Homo_sapiens.GRCh38.Ensembl99.genes.tsv", DataFrame, delim = "\t")
prot_coding = intersect(ensmbl[findall(ensmbl[:,"gene_biotype"] .== "protein_coding"), "SYMBOL"] , data[:,"Description"])
cds_data = data[in.(data[:,"Description"], (prot_coding,)), :]
values = Matrix(cds_data[:,3:end])
GTEX_tpm = Data(values, names(data)[3:end], cds_data[:, "Description"])
GTEX_tpm = Data(GTEX_tpm.data', GTEX_tpm.rows, GTEX_tpm.cols) # invert
# variances = vec(var(GTEX_tpm.data, dims = 1))
# keep = variances .> median(variances)
# GTEX_tpm_hv = Data(GTEX_tpm.data[:, keep], GTEX_tpm.cols[keep], GTEX_tpm.rows) 

# bson("data/gtex_tpm.bson", Dict("GTEX_tpm" => GTEX_tpm ))
#### TODO 
#### save to HDF5 !!!!
#### TODO 

pheno = CSV.read("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", DataFrame)
tissues = pheno[in.(pheno[:,"SAMPID"], (GTEX_tpm.rows,)), "SMTS"]
in.(pheno[:,"SAMPID"], (GTEX_tpm.rows,))
pheno
GTEX_tpm.rows
function Base.setindex!(f::HDF5.File, df::AbstractDataFrame, k::String) 
    g = create_group(f, k)
    for (name, vec) in pairs(eachcol(df))
        g[String(name)] = vec
    end
end
function DataFrames.DataFrame(g::HDF5.Group)
    convert(p) = (p.first, p.second[:]) # To pull data from the HDF5 dataset
    return DataFrame(Dict(map(convert, pairs(g))))
end
struct Bim
    h::String31
end 

f = h5open("GTEX.hdf5", "w")
f["data"] = GTEX_tpm.data
f["rows"] = GTEX_tpm.rows
f["cols"] = GTEX_tpm.cols
f["tissues"] = Array(tissues)
close(f)
# f["cds_data"] = cds_data

inf = h5open("GTEX.hdf5", "r")
data = log10.(inf["data"][:,:] .+ 1) 
rows = inf["rows"][:]
cols = inf["cols"][:]
tissues = inf["tissues"][:]
close(inf)


function label_binarizer(labels::Array)
    lbls = unique(labels)
    n = length(labels)
    m = length(lbls)
    binarizer = Array{Bool, 2}(undef, (n, m))
    for s in 1:n
        binarizer[s,:] = lbls .== labels[s]
    end 
    return binarizer
end 

using Flux
target = label_binarizer(tissues)'
loss(X, Y, model) = Flux.loss.MSE(model(X), Y)
model = Dense(125, 30, identity)
sign = rand(1: length(cols), 125)
X = data[:, sign]'
opt = Flux.ADAM(1e-2)
for e in 1:1000
    ps = Flux.params(model)
    l = Flux.Losses.mse(model(X), target)
    gs = gradient(ps) do
        Flux.Losses.mse(model(X), target)
    end
    Flux.update!(opt, ps, gs)
    println(l)
end 
preds = model(X) .== maximum(model(X), dims = 1)
preds && target
acc = target' .& preds'
pct = sum(acc) / n
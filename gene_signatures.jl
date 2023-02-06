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
using Random
using Dates
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

f = h5open("GTEX.out", "w")
f["data"] = GTEX_tpm.data
f["rows"] = GTEX_tpm.rows
f["cols"] = GTEX_tpm.cols
f["tissues"] = Array(tissues)
close(f)
# f["cds_data"] = cds_data

inf = h5open("GTEX.out", "r")
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

function accuracy(model, X, Y)
    n = size(X)[2]
    preds = model(X) .== maximum(model(X), dims = 1)
    acc = Y .& preds
    pct = sum(acc) / n
    return pct
end 

function train_logreg(X, Y; nepochs = 1000)
    model = gpu(Dense(size(X)[1], size(Y)[1], identity))
    opt = Flux.ADAM(1e-2)
    for e in 1:nepochs
        ps = Flux.params(model)
        l = Flux.logitcrossentropy(model(X), Y)
        # l = Flux.Losses.mse(model(X), target)
        gs = gradient(ps) do
            #Flux.Losses.mse(model(X), target)
            Flux.logitcrossentropy(model(X), Y)
        end
        Flux.update!(opt, ps, gs)
        # println(accuracy(model, X, Y))
    end
    return model 
end 

function split_train_test(X::Matrix, targets; nfolds = 5)
    folds = Array{Dict, 1}(undef, nfolds)
    nsamples = size(X)[1]
    fold_size  = Int(floor(nsamples / nfolds))
    ids = collect(1:nsamples)
    shuffled_ids = shuffle(ids)
    for i in 1:nfolds 
        tst_ids = shuffled_ids[collect((i-1) * fold_size +1: min(nsamples, i * fold_size))]
        tr_ids = setdiff(ids, tst_ids)
        train_x = X[tr_ids,:]
        train_y = targets[tr_ids, :]
        test_x = X[tst_ids, :]
        test_y = targets[tst_ids, :]
        folds[i] = Dict("train_x"=> train_x, "train_y" =>train_y, "test_x"=> test_x, "test_y" => test_y )
    end
    return folds  
end 
function set_dirs()
    outpath  = "./RES/SIGNATURES" # our output directory
    outdir = "$(outpath)/signatures_$(now())"
    mkdir(outdir)
    return outpath, outdir
end

targets = label_binarizer(tissues)

# 1,2,3,5,10,15,20,25,30,40,50,75,100,200,500
lengths = [1,2,3,5,10,15,20,25,30,40,50,75,100,200,500]
repn = 10
length_accs = Array{Float64, 2}(undef, (length(lengths) * repn, 2))

#### TCGA
basepath = "/home/muninn/scratch/golem-rpool-backup/golem-rpool/GDC"

function fetch_paths(path)
    dirs = readdir(path)
    paths = []
    tissues = Array{String,1}()
    samples = []
    for elem in dirs
        try 
            subpaths = readdir("$path/$elem")
            for subelem in subpaths
                if isfile("$path/$elem/$subelem/abundance.tsv")
                    push!(paths, "$path/$elem/$subelem/abundance.tsv")                
                    push!(tissues, elem)
                    push!(samples, subelem[1:8])
                end 
            end 
        catch err
        end 
    end 
    return paths, tissues, samples
end 

TCGA_paths, TCGA_tissues, TCGA_samples = fetch_paths(basepath)
template = CSV.read(TCGA_paths[1], DataFrame)
cols = template.target_id
rows = TCGA_samples
TCGA_TPM = Array{Float32, 2}(undef, (length(TCGA_paths), ngenes - 1))
using ProgressBars
for i::Int in ProgressBar(1:length(TCGA_samples)) 
    TCGA_TPM[i,:] = CSV.read(TCGA_paths[i], DataFrame).tpm
end 
# save in Data structure
TCGA_data = Data(TCGA_TPM, cols, rows)
f = h5open("TCGA.out", "w")
f["data"] = TCGA_data.data
f["rows"] = TCGA_data.rows
f["cols"] = TCGA_data.cols
f["tissues"] = Array(TCGA_tissues)
close(f)
# f["cds_data"] = cds_data

inf = h5open("TCGA.out", "r")
data = log10.(inf["data"][:,:] .+ 1) 
rows = inf["rows"][:]
cols = inf["cols"][:]
tissues = inf["tissues"][:]
close(inf)

ensmbl[findall(ensmbl[:,"gene_biotype"] .== "protein_coding"),:]
ENST = Array{String, 1}() 
raw = split.(cols,".")
[push!(ENST, raw[i][1]) for i in 1:size(raw)[1]]
ENST
TCGA_data = Data(data, ENST, rows)
gencode
variances = vec(var(data, dims = 1))
keep = variances .> sort(variances)[Int(round(length(variances) * 0.95))]
TCGA_TPM_hv = Data(data[:,keep], cols[keep], rows)
data = TCGA_TPM_hv.data
cols = TCGA_TPM_hv.cols
rows = TCGA_TPM_hv.rows
# variances = vec(var(GTEX_tpm.data, dims = 1))
# keep = variances .> median(variances)
# GTEX_tpm_hv = Data(GTEX_tpm.data[:, keep], GTEX_tpm.cols[keep], GTEX_tpm.rows) 


for (row, l) in enumerate(lengths)     
    for repl in 1:repn
        #loss(X, Y, model) = Flux.loss.MSE(model(X), Y)
        sign = rand(1: length(cols), l)
        X = data[:, sign]
        folds = split_train_test(X, targets)
        accs = []
        for (foldn, fold) in enumerate(folds)
            train_x = gpu(fold["train_x"]')
            train_y = gpu(fold["train_y"]')
            test_x = gpu(fold["test_x"]')
            test_y = gpu(fold["test_y"]')

            model = train_logreg(train_x, train_y, nepochs = 1000)
            println("Length $l Rep $repl Fold $foldn Train : ", accuracy(model, train_x, train_y))
            println("Length $l Rep $repl Fold $foldn Test  : ", accuracy(model, test_x, test_y))
            push!(accs, accuracy(model, test_x, test_y))
        end
        length_accs[(row - 1) * repn + repl,:] =  Array{Float64}([l, mean(accs)])
    end 
end 
df = DataFrame(Dict([("lengths", length_accs[:,1]), ("tst_acc", length_accs[:,2])]))
outpath, outdir = set_dirs()
CSV.write("$outdir/TCGA_tst_accs.csv", df)

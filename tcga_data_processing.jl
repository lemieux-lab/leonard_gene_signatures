# load manifests 
using CSV
using DataFrames
using JSON
using ProgressBars
using HDF5
using Statistics
CLIN_FULL = CSV.read("data/tcga_clinical_raw.tsv", DataFrame)
MANIFEST = CSV.read("data/gdc_manifest_GE_2023-02-02.txt", DataFrame)
IDS = CSV.read("data/sample.tsv", DataFrame)
FILES = "data/files.2023-02-06.json"
baseurl = "https://api.gdc.cancer.gov/data"
basepath = "data/TCGA"
J = JSON.parsefile(FILES)
features = ["case_id", "project_id", "gender", "age_at_index","age_at_diagnosis", "days_to_death", "days_to_last_follow_up", "primary_diagnosis", "treatment_type"]
CLIN_FULL
CLIN = CLIN_FULL[:, features]

sum(CLIN.days_to_death .!= "'--") # UNCENSORED data
nfiles = length(J)
case_ids = [elem["cases"][1]["case_id"] for elem in J]
findall(case_ids .== "1e10278b-1367-4a3f-b2ad-89566bbc3e5c")

f = open("data/fetch_data.sh", "w")
## FECTHING DATA 
for i::Int in ProgressBar(2838:nfiles)
    file_id = J[i]["file_id"]
    println(file_id)
    case_id = J[i]["cases"][1]["case_id"]
    println(case_id)
    cmd = "curl $baseurl/$file_id -o $basepath/$case_id\n"
    write(f, cmd)
    # cmd = `curl $baseurl/$file_id -o $basepath/$case_id`
    #run(cmd)
end 
close(f)
## for filename / fileID 
## fetch file 
## TCGA/case_id/fname.tsv
### MERGING and STORING to HDF5


sample_data = CSV.read("$basepath/$(case_ids[1])", DataFrame, delim = "\t", header = 2)
sample_data = sample_data[5:end, ["gene_name", "tpm_unstranded"]]

struct TCGA_data
    data::Matrix
    rows::Array 
    cols::Array
end 
ngenes = size(merged)[1]
nsamples = length(case_ids)
m=Array{Float32, 2}(undef, (nsamples, ngenes))

for i in ProgressBar(1:5500)
    data = CSV.read("$basepath/$(case_ids[i])", DataFrame, delim = "\t", header = 2)
    data = data[5:end, ["gene_name", "tpm_unstranded"]]
    m[i, :] = data.tpm_unstranded
end
tmp_data = TCGA_data(m, case_ids, Array{String}(sample_data.gene_name))
# HDF5
# writing to hdf5 
f = h5open("data/TCGA_STAR_TPM.h5", "w")
f["data"] = tmp_data.data
f["rows"] = tmp_data.rows
f["cols"] = tmp_data.cols
close(f)

#loading h5 file 
inf = h5open("data/TCGA_STAR_TPM.h5", "r")
tpm_data = log10.(inf["data"][:,:] .+ 1)
case_ids = inf["rows"][:]
gene_names = inf["cols"][:] 
close(inf)
tpm_data
# subset
keep1 = findall(vec(sum(tpm_data, dims=2 )) .!= 0)
cases = case_ids[keep1]
tcga1 = TCGA_data(tpm_data[keep1,:], cases, gene_names)

# intersect with clin data 
uniq_case_id = unique(CLIN.case_id)
keep2 = [in(c,uniq_case_id ) for c in cases]
tcga2 = TCGA_data(tcga1.data[keep2,:], tcga1.rows[keep2], tcga1.cols)

# map to tissues 
cid_pid = Dict([(cid, pid) for (cid, pid) in zip(CLIN.case_id, Array{String}(CLIN.project_id))])
tissues = [cid_pid[c] for c in tcga2.rows]
vars = vec(var(tcga.data, dims = 1))  
hv = vec(var(tcga.data, dims =1 )) .> sort(vars)[Int(round(0.75 * ngenes))]

tcga_hv = TCGA_data(tcga2.data[:,hv], tcga2.rows, tcga2.cols[hv])

f = h5open("data/TCGA_STAR_TPM_subset.h5", "w")
f["data"] = tcga_hv.data
f["rows"] = tcga_hv.rows
f["cols"] = tcga_hv.cols
f["tissues"] = tissues
close(f)

using TSne
@time TCGA_tsne = tsne(tcga_hv.data, 2, 50, 1000, 30.0;verbose=true,progress=true)

using CairoMakie
using AlgebraOfGraphics

TSNE_df = DataFrame(Dict("dim_1" => TCGA_tsne[:,1], "dim_2" => TCGA_tsne[:,2], "tissue" => tissues))

q = AlgebraOfGraphics.data(TSNE_df) * mapping(:dim_1, :dim_2, color = :tissue, marker = :tissue) * visual(markersize = 15,strokewidth = 0.5, strokecolor =:black)

main_fig = draw(q ; axis=(width=1024, height=1024,
                title = "2D TSNE by tissue type on TCGA data, number of input genes: $(size(tcga_hv.data)[2]), nb. samples: $(size(tcga_hv.data)[1])",
                xlabel = "TSNE 1",
                ylabel = "TSNE 2"))
save("RES/TCGA_5500_samples.png", main_fig, pt_per_unit = 2)


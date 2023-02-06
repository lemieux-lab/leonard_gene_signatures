# load manifests 
using CSV
using DataFrames
using JSON
CLIN = CSV.read("data/tcga_clinical_raw.tsv", DataFrame)
MANIFEST = CSV.read("data/gdc_manifest_GE_2023-02-02.txt", DataFrame)
IDS = CSV.read("data/sample.tsv", DataFrame)
FILES = "data/files.2023-02-06.json"
baseurl = "https://api.gdc.cancer.gov/data"
basepath = "data/TCGA"
J = JSON.parsefile(FILES)
J[1]["cases"][1]["case_id"]
features = ["gender", "age_at_index","age_at_diagnosis", "days_to_death", "days_to_last_follow_up", "primary_diagnosis", "treatment_type"]
CLIN = CLIN[:, features]

sum(CLIN.days_to_death .!= "'--") # UNCENSORED data

## FECTHING DATA 
for i in 1:10
    file_id = J[i]["file_id"]
    println(file_id)
    case_id = J[i]["cases"][1]["case_id"]
    println(case_id)
    cmd = `curl $baseurl/$file_id -o $basepath/$case_id`
    run(cmd)
end 
## for filename / fileID 
## fetch file 
## TCGA/case_id/fname.tsv
### MERGING and STORING to HDF5
case_ids = [elem["cases"][1]["case_id"] for elem in J]


merged = CSV.read("$basepath/$(case_ids[1])", DataFrame, delim = "\t", header = 2)
merged = merged[5:end, ["gene_name", "tpm_unstranded"]]
rename!(merged, ["gene_name", "$(case_ids[1])"])
merged
for case_id in case_ids[2:10]
    data = CSV.read("$basepath/$case_id", DataFrame, delim = "\t", header = 2)
    data = data[5:end, ["gene_name", "tpm_unstranded"]]
    rename!(data, ["gene_name", "$case_id"])
    sum(merged.gene_name .!= data.gene_name)
    #merged = leftjoin(merged, data, on = :gene_name)
end 
merged
mkdir("data/GDC")
MANIFEST.id
MANIFEST.id[1]
l = [in(elem, MANIFEST.id) for elem in unique(IDS.sample_id)]
sum(l)
CLIN[CLIN.case_id .== MANIFEST.id[1],:]
for i in 1:10
    id = MANIFEST.id[i]
    cmd = `curl $baseurl/$id > data/GDC/$id`
    println(cmd)
end
unique(CLIN.case_id)
MANIFEST.id
tmp = CLIN[CLIN.]

# hdf5 ? 
using DataFrames, CSV, Tables

d = "/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Padded Data"

files = readdir(d, join=true)

filter!(f -> occursin("Reference", f), files)


sort!(files)

A = []
for f in files
        data = open(f) do fh
                parse.(Float64, readlines(fh))
        end
        append!(A, data)
end
B = reshape(A, 250, :)



CSV.write("/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/Reference.csv", Tables.table(B))

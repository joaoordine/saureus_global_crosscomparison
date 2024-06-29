# Pan-resistome inferral ------- PanViTa v1.0 - clustermap with data grouped with Euclidean distance measure  

## Installing tool and its dependencies 
git clone https://github.com/dlnrodrigues/panvita.git
cd panvita
python3 panvita.py -u
python3 panvita.py -h


## Moving all gbk files from Prokka into a single directory 
cd prokka_output
mkdir all_gbk_files
for dir in GCA_*; do
    scp /temporario2/11217468/projects/saureus_global/prokka_output/"$dir"/*.gbk all_gbk_files && 
    echo "Copy of $dir was sent"
done 

## If any dependencie isn't installed automatically, run the following
###conda install -c anaconda wget
###conda install seaborn
###conda install pandas
###conda install -c conda-forge matplotlib
###conda install -c anaconda basemap

## Run the tool 
input_dir="/temporario2/11217468/projects/saureus_global/prokka_output/all_gbk_files"

python3 panvita.py -card -vfdb -bacmet -i 80 -c 80 "$input_dir"/*.gbk 





mkdir B6_cellphoneDB
cd B6_cellphoneDB
cellphonedb method statistical_analysis ../B6_filtered_meta.txt ../B6_filtered_hcount.txt 
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot ../B6_filtered_meta.txt

cd ..
mkdir SLE.yaa_cellphoneDB
cd SLE.yaa_cellphoneDB
cellphonedb method statistical_analysis ../SLE.yaa_filtered_meta.txt ../SLE.yaa_filtered_hcount.txt 
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot ../SLE.yaa_filtered_meta.txt

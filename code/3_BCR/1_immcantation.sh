
docker run -it --workdir /data -v "$(pwd):/data:z" immcantation/suite:4.1.0 bash
changeo-10x -s filtered_contig.fasta -a filtered_contig_annotations.csv -o . \
    -g mouse -t ig -x auto
light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e filtered_contig_light_productive-T.tsv \
    -o filtered_contig_heavy_clone-light.tsv
CreateGermlines.py -d filtered_contig_heavy_clone-light.tsv -g dmask --cloned \
    -r /usr/local/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta \
    /usr/local/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
    /usr/local/share/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
    --outname filtered_contig_heavy
BuildTrees.py -d filtered_contig_heavy_germ-pass.tsv --minseq 3 --clean all \
    --igphyml --collapse --nproc 2 --asr 0.9
exit






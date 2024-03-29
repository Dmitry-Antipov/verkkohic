#!/bin/bash
#module load gcc/9.2.0 
#module load python/3.8
#module load mashmap
#module load snakemake
#PSTOOLS=/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/
#VERKKO=/data/korens/devel/verkko-tip/

HIC1=$3/*hic*/*R1*fastq.gz
HIC2=$3/*hic*/*R2*fastq.gz
#verkko_output_folder script_output hic-reads
display_usage() { 
    echo -e "\nUsage: $0 [previous verkko run] [output folder] [reads folder] \n" 
    } 
# if less than two arguments supplied, display usage 
if [  $# -le 1 ] 
   then 
   display_usage
   exit 1
fi 



echo "---Running consensus on graph edges to get homopolymer uncompressed seqs"
sh $VERKKO/bin/verkko  --paths $1/6-layoutContigs/consensus_paths.txt --assembly $1 -d $2/consensus_unitigs/ --hifi $3/hifi/*fast*.gz --nano $3/ont/*fast*.gz

echo "---Preprocessing graph files"
mkdir -p $2
#THIS DOES NOT WORK WITH TRIO RUNS(
#cat $1/6-layoutContigs/unitig-popped.layout.scfmap | awk '{if (match($1, "path")) print $2"\t"$3}' > $2/contigs_rename.map 
#python3 $VERKKO/lib/verkko/scripts/process_reads.py rename $2/assembly.fasta $2/contigs_rename.map $1/assembly.fasta

#cp $2/consensus_unitigs/assembly.fasta $2/unitigs.fasta

cat $2/consensus_unitigs/6-layoutContigs/unitig-popped.layout.scfmap | awk '{if (match($1, "path")) print $2"\t"$3}' > $2/contigs_rename.map 
python3 $VERKKO/lib/verkko/scripts/process_reads.py rename $2/unitigs.fasta  $2/contigs_rename.map $2/consensus_unitigs/assembly.fasta

awk '/^S/{print ">"$2"\n"$3}' $1/assembly.homopolymer-compressed.gfa | fold > $2/unitigs.hpc.fasta
cp $1/assembly.homopolymer-compressed.noseq.gfa $2/unitigs.hpc.noseq.gfa


echo "---Running msahmap"
#homopolymer compressed unitigs for mashmap
mashmap -r $2/unitigs.hpc.fasta -q $2/unitigs.hpc.fasta -t 8 -f none --pi 95 -s 10000 -o $2/mashmap.out
cat $2/mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > $2/unitigs.matches

echo "---Mapping reads with pstools"
#hicmapping
echo "Mapping $HIC1 and $HIC2 to $2/unitigs.fasta"
echo "pstools will write BIG temporary files to current directory $PWD"
if [ ! -e $2/map_uncompressed.out ]; then
   echo "$PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/unitigs.fasta <(zcat $HIC1) <(zcat $HIC2)"
   $PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/unitigs.fasta <(zcat $HIC1) <(zcat $HIC2)
   mv hic_name_connection.output $2/hic_mapping.byread.output
fi


if [ -n "$SLURM_JOB_ID" ] ; then
THEPATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
echo "THEPATH1 $THEPATH"
THEPATH=`echo "${THEPATH}" | head -1`
echo "THEPATH2 $THEPATH"

echo "Detected SLURM"
else
THEPATH=$(realpath $0)
fi

echo "Path of script detected $THEPATH"
SCRIPT_DIR=$(dirname $THEPATH)
SCRIPT_DIR=$(realpath $SCRIPT_DIR)
SCRIPT_DIR=`echo "${SCRIPT_DIR}" | head -1`
echo "---Running main script"
echo "Directory of the script $SCRIPT_DIR "
echo "Path to script..."
echo "$SCRIPT_DIR/hicverkko.py"

python3 $SCRIPT_DIR/hicverkko.py $1 $2 


echo "---Running rukki on the resulting clustering"
echo "---Running rukki on the resulting clustering"
params=""
params="$params --init-assign $2/out_init_ann.csv"
params="$params --refined-assign $2/out_refined_ann.csv"
params="$params --final-assign $2/out_final_ann.csv"
params="$params --marker-sparsity 5000"
params="$params --issue-sparsity 1000"
params="$params --try-fill-bubbles"
params="$params --fillable-bubble-diff 1000"
params="$params --fillable-bubble-len 500000"
params="$params --assign-tangles --tangle-allow-deadend"
params="$params --issue-ratio 1."
params="$params --solid-homozygous-cov-coeff 1.1"
params="$params --solid-ratio 1.5"
params="$params --hap-names haplotype1,haplotype2"

if [ xtrio = xtrio ]; then
   params="$params --marker-ratio 5."
else
   params="$params --marker-ratio 3."
fi

$VERKKO/lib/verkko/bin/rukki trio -g $2/unitigs.hpc.noseq.gfa -m $2/hicverkko.colors.tsv              -p $2/rukki.paths.tsv $params
$VERKKO/lib/verkko/bin/rukki trio -g $2/unitigs.hpc.noseq.gfa -m $2/hicverkko.colors.tsv --gaf-format -p $2/rukki.paths.gaf $params

echo "---Running final verkko consensus on paths"
sh $VERKKO/bin/verkko  --paths $2/rukki.paths.gaf --assembly $1 -d $2/final_consensus/ --hifi $3/hifi/*fast*.gz --nano $3/ont/*fast*.gz

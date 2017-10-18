#1/bin/bash

   sra_dir=$1
   t=$2
   fastq-dump --split-files ${sra_dir}"/"*.sra -O ${sra_dir}
   cat *_1.fastq > "${basename}-1.fastq" &
   cat *_2.fastq > "${basename}-2.fastq"
   rm *_?.fastq
   fadir="/home/GC/UCSC/Sequences/chromFa.mm10"   fa="${fadir}/chr[1-9,M,X,Y].fa ${fadir}/chr1[0-9]"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++"
echo ""
/home/GC/software/segemehl/segemehl.x -i /home/GC/index/segemehl-index/BS-paired/mm10/chr.ctidx -j /home/GC/index/segemehl-index/BS-paired/mm10/chr.gaidx -d ${fa} -q /home/lysong/sperm_MethylC-Seq/SRR1286765-73_1.fastq -p /home/lysong/sperm_MethylC-Seq/SRR1286765-73_1.fastq -F 1 -M 1 -H 1 -t ${t} -o /home/lysong/sperm_MethylC-Seq/sperm.sam  


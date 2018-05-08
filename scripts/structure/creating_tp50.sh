~/share/Computing/htslib-master/tabix -fh ~/share/ADMIXTURE/filteredNew.vcf.bgz Scaffold0:0-1 | bgzip > ~/share/ADMIXTURE/neutral/neutral.vcf.bgz
xargs -a ~/share/ADMIXTURE/neutral/neutral.txt -I {} ~/share/Computing/htslib-master/tabix -f ~/share/ADMIXTURE/filteredNew.vcf.bgz {} | bgzip >> ~/share/ADMIXTURE/neutral/neutral.vcf.bgz

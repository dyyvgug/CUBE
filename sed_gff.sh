# deal with GFF files, to extract cytosolic RP genes.
sed -i '/ribosomal/!d' *.gff
sed -i '/apicoplast/d' *.gff
sed -i '/mitochondrial/d' *.gff
sed -i '/ubiquitin/d' *.gff
sed -i '/CDS/!d' *.gff

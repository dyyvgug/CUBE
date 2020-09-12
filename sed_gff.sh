sed -i '/ribosomal/!d' *.gff
sed -i '/apicoplast/d' *.gff
sed -i '/mitochondrial/d' *.gff
sed -i '/ubiquitin/d' *.gff
sed -i '/CDS/!d' *.gff

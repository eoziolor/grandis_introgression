cat top50_annot | perl -lne 'print $1 if /gene=\s*(.*)/' | grep -v "uncharacterized" | sed "s/;.*//" | grep -v "LOC" | uniq

cat top50_annot | perl -lne 'print $1 if /gene=\s*(.*)/' | grep -v "uncharacterized" | grep "LOC" | perl -lne 'print $1 if /product=\s*(.*)/' | sed "s/;.*//" | uniq

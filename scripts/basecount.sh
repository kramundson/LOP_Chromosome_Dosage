for read in `ls *.fq` 
do
	echo -n "$read " ; cat $read | paste - - - - | cut -f 2 | wc -c
done > basecount.txt

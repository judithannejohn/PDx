echo "Number of alignments of reads to reference sequences in the file"
cut -f10 sequence.txt.sam | sort | uniq | wc -l
echo "Ans:99952"
echo "99952/100000 = 0.99% reads participating in the alignment"

echo "uniquely mapped"
grep -c XT:A:U sequence.txt.sam
echo "Ans: 94035"

echo "multi hit reads"
grep -c XT:A:R sequence.txt.sam
echo "Ans: 5518"

echo "Number of reads wrongly aligned"
echo "(No. of participated reads for alignment)- (uniqly mapped + multi hit map)"
echo "Ans: 99952-(94035-5518)= 399"
echo "error rate= 399/99952"
echo "0.39% wrongly mapped"

PEAKS=(
    "../2018-11-27_peak-filtering/Filter/logq_5.0/Ctrl_Rep1.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/Ctrl_Rep2.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/Ctrl_Rep3.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/MB6_Rep1.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/MB6_Rep2.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/MB6_Rep3.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/MB6Pen_Rep1.bed"
    "../2018-11-27_peak-filtering/Filter/logq_5.0/MB6Pen_Rep2.bed"
)

echo "Condition-based consensus peaks"
echo "  Ctrl"
cat ${PEAKS[@]:0:3} | sort -k1,1 -k2,2n | bedtools merge -i stdin > Consensus/Ctrl.bed
bedtools intersect -a Consensus/Ctrl.bed -b ${PEAKS[@]:0:3} -c -sorted > Consensus/Ctrl.counts.bed
echo "  MB6"
cat ${PEAKS[@]:3:3} | sort -k1,1 -k2,2n | bedtools merge -i stdin > Consensus/MB6.bed
bedtools intersect -a Consensus/MB6.bed -b ${PEAKS[@]:3:3} -c -sorted > Consensus/MB6.counts.bed
echo "  MB6Pen"
cat ${PEAKS[@]:6:2} | sort -k1,1 -k2,2n | bedtools merge -i stdin > Consensus/MB6Pen.bed
bedtools intersect -a Consensus/MB6Pen.bed -b ${PEAKS[@]:6:2} -c -sorted > Consensus/MB6Pen.counts.bed

# keep sites that exist in at least 2/3 replicates
echo "Filtering sites"
awk '{if ($4 > 1) print $0}' Consensus/Ctrl.counts.bed > Consensus/Ctrl.filtered.bed
awk '{if ($4 > 1) print $0}' Consensus/MB6.counts.bed > Consensus/MB6.filtered.bed
awk '{if ($4 > 1) print $0}' Consensus/MB6Pen.counts.bed > Consensus/MB6Pen.filtered.bed

# generate consensus for all 3 conditions from the intersections
echo "Generating 3-condition consensus"
cat Consensus/{Ctrl,MB6,MB6Pen}.filtered.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > Consensus/consensus.all.bed

# generate consensus for all 3 conditions from the intersections
echo "Generating 2-condition consensuses"
cat Consensus/{Ctrl,MB6}.filtered.bed | sort -k1,1 -V -k2,2n | bedtools merge -i stdin > Consensus/consensus.Ctrl-MB6.bed
cat Consensus/{Ctrl,MB6Pen}.filtered.bed | sort -k1,1 -V -k2,2n | bedtools merge -i stdin > Consensus/consensus.Ctrl-MB6Pen.bed
cat Consensus/{MB6,MB6Pen}.filtered.bed | sort -k1,1 -V -k2,2n | bedtools merge -i stdin > Consensus/consensus.MB6-MB6Pen.bed

echo "Done"

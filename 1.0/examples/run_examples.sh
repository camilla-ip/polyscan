#!/usr/bin/env bash

# run_examples.sh

# Shared input files

REFNAME=A_X02763
REFPROT=./input/references/${REFNAME}_prot.fasta
REFNUCL=./input/references/${REFNAME}_nucl.fasta
SAMPLETABLE=./input/sample-list.tsv
AAALNDIR=./input/protein_alignments
NTALNDIR=./input/nt_alignments

# Run examples 1 to 2

echo "Starting example runs"
echo
echo '-----'
for i in {1..2} ; do

    POLYTABLE=./input/polymorphism-list-example${i}.tsv
    OUTDIR=./output/Ex${i}
    OUTPREFIX=Ex${i}
    RUNOUT=${OUTDIR}/${OUTPREFIX}.run.summary
    RUNDIFF=${OUTDIR}/${OUTPREFIX}.run.diff

    echo
    echo "Running Example ${i}"
    set -x
    ../polyscan.py ${REFNAME} ${REFPROT} ${REFNUCL} ${POLYTABLE} ${SAMPLETABLE} ${AAALNDIR} ${NTALNDIR} ${OUTDIR} ${OUTPREFIX} &> ${RUNOUT}

    if [ -f ${RUNDIFF} ] ; then rm ${RUNDIFF} ; fi
    touch ${RUNDIFF}
    diff ${OUTDIR}/expected/${OUTPREFIX}_sample_polymorphisms.tsv ${OUTDIR}/${OUTPREFIX}_sample_polymorphisms.tsv &> ${RUNDIFF}
    diff ${OUTDIR}/expected/${OUTPREFIX}_summary_polymorphisms.tsv ${OUTDIR}/${OUTPREFIX}_summary_polymorphisms.tsv &> ${RUNDIFF}
    diff ${OUTDIR}/expected/${OUTPREFIX}_summary_polymorphisms_pvalue.tsv ${OUTDIR}/${OUTPREFIX}_summary_polymorphisms_pvalue.tsv &> ${RUNDIFF}
    set +x
    echo "See ${RUNOUT} for any run-time warnings or errors."
    echo "See ${RUNDIFF} for any differences between the expected and test output. If the file is empty, there are no errors."
    echo
    echo '-----'

done
echo "Finished example runs"

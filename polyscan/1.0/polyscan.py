#!/usr/bin/env python3

from Bio import SeqIO
import os
import pandas as pd
from scipy.stats import fisher_exact
import sys

# Global variables

_version = '1.0'

_usage_msg = '''
Usage:   {prog} [9 arguments]
Version: {ver}

For each sample, scan for the specified polymorphisms, and output a table of polymorphisms in each sample.

Tabulate the number of samples, grouped by HCC or non-HCC diagnosis, that contain each
polymorphism, and the p-value for the Fisher's Exact Test for the HCC vs non-HCC counts.

Arguments
---------
            refname : The reference genome name used in input file IDs (e.g., A_X02763)
  refprot_fastafile : Amino-acid sequences with id=REFNAME_PROTNAME [FASTA]
  refnucl_fastafile : Nucleotide genome with id=REFNAME [FASTA]
  polytable_tsvfile : Table of polymorphisms with 1 header row and fields
                      poly_ prefix then [_group, _name, _desc, _level, _gene, _type, _refval _refpos, _altval]
                      Multiple rows for the same (poly_group, poly_name) means the criteria specified in all rows
                      must be present for the (poly_group, poly_name) variant to be considered present[TSV]
sampletable_tsvfile : Table of samples with 1 header row and fields
                      [dataid, samplename, studynumber, diagnosis, HCC] [TSV]
          aaaln_dir : Directory containing DATAID_PROTNAME_AA_align_REFNAME.fasta files
                      Each file contains an alignment of the reference protein (>REFNAME_PROTNAME) and
                      the sample protein (>DATAID_bam_ref_PROTNAME).
                      If a file for a particular PROTNAME does not exist, the protein could not be recovered
                      from the de novo assembly for that sample.
          ntaln_dir : Directory containing DATAID_nt_aligned.fasta files
                      Each file contains an alignment of the reference genome (>REFNAME) and
                      the de novo assembly for the sample (>DATAID)
        outdata_dir : Output directory
     outdata_prefix : Prefix for output files

'''.format(prog=os.path.basename(sys.argv[0]), ver=_version)

_refname = None
_refprot_fastafile = None
_refnucl_fastafile = None
_polytable_tsvfile = None
_sampletable_tsvfile = None
_aaaln_dir = None
_ntaln_dir = None
_outdata_dir = None
_outdata_prefix = None

_debug = False
_num_dec_pl = 1
_Fishers_Exact_Test_p_value_signif_threshold = 0.05

# Data structures

'''
_refprot_fastafile

    A[protname] = seq
    protnameL = [protname, ...]

_refnucl_fastafile

    N[refname] = seq

_polytable_tsvfile

    P[(poly_group, poly_name)] = {
        'desc': X,
        'spec': [ { 'level': X, 'gene': X, 'type': X, 'refval': [X,...], 'refpos': [X,...], 'altval': [X,...] }, ... ]
    }
    polygroupL = [poly_group , ...]
    polykeyL = [(poly_group, poly_name), ...]

_sampletable_tsvfile

    S[dataid] = { "samplename": X, "studynumber": X, "diagnosis": X, "hcc": X }
    sampleL = [dataid, ...]

_aaaln_dir

    I[(dataid, protname)] = { 'ref': seq, 'sample': seq }

_ntaln_dir

    J[dataid] = { 'ref': seq, 'sample': seq }

'''

# Parse input files

def Parse_Fasta(fastafile, unwantedprefix=""):
    'Return the contents of the FASTA file as D[id] = seq and the list of ids in L.'
    D = {}
    L = []
    for record in SeqIO.parse(fastafile, 'fasta'):
        id = record.id[len(unwantedprefix):] \
            if len(unwantedprefix) and record.id.startswith(unwantedprefix) \
            else record.id
        D[id] = str(record.seq).upper()
        L.append(id)
    return D, L

def Parse_PolyTable(polytablefile):
    '''
    Return the contents of the polytable as
    P[(group, name)]['desc'] = X
    P[(group, name)]['spec'] = [ { 'level': X, 'gene': X, 'type': X, 'refval': X, 'refpos': X, 'altval': X }, ... ]
    polygroupL = [group, ...]
    polykeyL = [(group, name), ...]
    '''
    P = {}
    polygroupL = []
    polykeyL = []
    with open(polytablefile, 'r') as in_fp:
        linecnt = 0
        for line in in_fp:
            linecnt += 1
            if linecnt == 1:
                continue
            L = line.rstrip('\n').rstrip('\r').split('\t')
            poly_group, poly_name, poly_desc, poly_level, poly_gene, poly_type, poly_refval, poly_refpos, poly_altval = L
            if not poly_group in P:
                P[poly_group] = {}
                polygroupL.append(poly_group)
            key = (poly_group, poly_name)
            if not key in polykeyL:
                polykeyL.append(key)
            if not poly_name in P[poly_group].keys():
                P[poly_group][poly_name] = { 'desc' : poly_desc, 'spec': [] }
            if poly_type == 'SNP':
                P[poly_group][poly_name]['spec'].append({
                    'level': poly_level, 'gene': poly_gene, 'type': poly_type,
                    'refval': poly_refval.split(','),
                    'refpos': [int(x) for x in poly_refpos.split('-')],
                    'altval': poly_altval.split(',') })
            elif poly_type == 'DELat':
                P[poly_group][poly_name]['spec'].append({
                    'level': poly_level, 'gene': poly_gene, 'type': poly_type,
                    'refval': None, 'refpos': [int(x) for x in poly_refpos.split('-')], 'altval': None})
            elif poly_type == 'DELany':
                P[poly_group][poly_name]['spec'].append({
                    'level': poly_level, 'gene': poly_gene, 'type': poly_type,
                    'refval': None, 'refpos': None, 'altval': None})
    return P, polygroupL, polykeyL

def Parse_SampleTable(samplefile):
    'Return S[dataid] = { "samplename": X, "studynumber": X, "diagnosis": X, "hcc": X } and sampleL = [dataid, ...].'
    S = {}
    sampleL = []
    with open(samplefile, 'r') as in_fp:
        linecnt = 0
        for line in in_fp:
            linecnt += 1
            if linecnt == 1:
                continue
            L = line.rstrip('\n').split('\t')
            dataid, samplename, study_number, diagnosis, hcc = L
            hcc = int(hcc)
            S[dataid] = { 'samplename': samplename, 'studynumber': study_number, 'diagnosis': diagnosis, 'hcc': hcc }
            sampleL.append(dataid)
    return S, sampleL

def Sample_Prot_Aln_Path(indir, dataid, protname, refname):
    return os.path.join(indir, dataid+'_'+protname+'_AA_align_'+refname+'.fasta')

def Sample_Nt_Aln_Path(indir, dataid):
    return os.path.join(indir, dataid+'_nt_aligned.fasta')

def Parse_ProtAlignmentFiles(S, dataidL, protnameL, indir, refname):
    '''
    Return the alignments of ref vs sample proteins in I[(dataid, protname)] = { "ref": seq, "sample": seq }
    If there is no alignment for that gene in that sample, I[(dataid, protname)] = None
    '''
    I = {}
    for dataid in dataidL:
        for protname in protnameL:
            key = (dataid, protname)
            aln_path = Sample_Prot_Aln_Path(indir, dataid, protname, refname)
            if os.path.exists(aln_path) and os.path.getsize(aln_path) > 0:
                D, L = Parse_Fasta(aln_path)
                refaln = D[refname+'_'+protname]
                samplealn = D[dataid+'_bam_ref_'+protname]
                I[key] = { 'ref': refaln, 'sample': samplealn }
            else:
                I[key] = None
    return I

def Parse_NtAlignmentFiles(S, indir, refname):
    'Store contents of input files with >REFNAME and >DATAID records, in J[dataid] = { "ref": seq, "sample": seq }'
    J = {}
    dataidL = list(S.keys())
    dataidL.sort()
    for dataid in dataidL:
        aln_path = Sample_Nt_Aln_Path(indir, dataid)
        if os.path.exists(aln_path) and os.path.getsize(aln_path) > 0:
            F, L  = Parse_Fasta(aln_path)
            J[dataid] = { 'ref': F[refname], 'sample': F[dataid] }
        else:
            J[idataid] = None
    return J

def Polytable_Row_Ok(poly_group, poly_name, E, refname, A, N):
    '''
    Return 1 if the polymorphism spec makes sense. Rows of type:
       SNP : must have a csv string of refvals
     DELat : must have a refval of '-' and either a single number N or a range N-N for refpos
    DELany : must have a refval of '-' and refpos of '-'
    If the reference nt/aa is not in the list of refpos values, print a warning.
    '''

    if E['type'] == 'SNP':
        refvalL = E['refval']
        if E['level'] == 'aa':
            refval = A[E['gene']][E['refpos'][0]-1]
        else:
            refval = N[refname][E['refpos'][0]-1]
        if refval not in refvalL:
            error_msg = 'Warning: Polytable SNP refval disagrees with refprot FASTA ({g}, {n}): table={a} fasta={b}'.format(
                g=poly_group, n=poly_name, a=','.join(E['refval']), b=refval)
            return False, error_msg
        return True, None
    elif E['type'] == 'DELat':
        return True, None
    elif E['type'] == 'DELany':
        return True, None
    else:
        error_msg = 'Error: Invalid poly_type detected *{t}*'.format(t=E['type'])
        return False, error_msg

def Validate_Polytable_vs_ReferenceSequences(A, P, polygroupL, polykeyL, N, refname):
    '''
    Check that the (refval+refpos) values in the polytable agree with the amino-acid sequence values
    in the reference FASTA files.
    '''
    for (poly_group, poly_name) in polykeyL:
        if _debug:
            sys.stdout.write('Debug: Validate_Polytable_vs_ReferenceSequences: poly_group={g}, poly_name={n}\n'.format(g=poly_group, n=poly_name))
        for E in P[poly_group][poly_name]['spec']:
            ok, error_msg =  Polytable_Row_Ok(poly_group, poly_name, E, refname, A, N)
            if not ok:
                sys.stdout.write(error_msg+'\n')

# Compute Values

def AlnPos2Idx(alnS, pos):
    'Convert 1-based position to 0-based index of a string with gaps.'
    nogapsS = alnS.replace('-', '')
    if pos < 1 or pos > len(nogapsS):
        return -1
    idx = 0
    cnt = 1
    while idx < len(alnS) and cnt < pos:
        idx += 1
        cnt += (alnS[idx] != '-')
    return idx

def Sample_SNP_Result(dataid, poly_group, poly_name, E, I, J):
    if E['level'] == 'aa':
        if not (dataid, E['gene']) in I or I[(dataid, E['gene'])] is None:
            return 'NA'
        ref_aln_seq = I[(dataid, E['gene'])]['ref']
        sample_aln_seq = I[(dataid, E['gene'])]['sample']
    else:
        if not dataid in J or J[dataid] is None:
            return 'NA'
        ref_aln_seq = J[dataid]['ref']
        sample_aln_seq = J[dataid]['sample']

    aln_idx = AlnPos2Idx(ref_aln_seq, E['refpos'][0])
    sample_variant = sample_aln_seq[aln_idx]
    if '*' in E['altval']:
        result = int((aln_idx+1) > len(sample_aln_seq))
    else:
        result = int(sample_variant in E['altval'])
    if _debug:
        if result == 1:
            sys.stdout.write('Debug: Sample_SNP_Result {d}, {g}, {n}, E={e}, sample_variant={v} => {r}\n'.format(
                d=dataid, g=poly_group, n=poly_name, e=E, v=sample_variant, r=result))
    return result

def Sample_DELat_Result(dataid, poly_group, poly_name, E, I, J):
    'Return 1 if the aligned sample sequence has a deletion at the specified position range.'

    if E['level'] == 'aa':
        if not (dataid, E['gene']) in I or I[(dataid, E['gene'])] is None:
            return 'NA'
        ref_aln_seq = I[(dataid, E['gene'])]['ref']
        sample_aln_seq = I[(dataid, E['gene'])]['sample']
    else:
        if not dataid in J or J[dataid] is None:
            return 'NA'
        ref_aln_seq = J[dataid]['ref']
        sample_aln_seq = J[dataid]['sample']

    startidx = AlnPos2Idx(ref_aln_seq, E['refpos'][0])
    endidx = AlnPos2Idx(ref_aln_seq, E['refpos'][1])
    if (endidx >= len(ref_aln_seq)):
        sys.stderr.write('Error: Sample_DELat_Result: {d}, {g}, {n}, E={e} - endidx after end of ref_aln_seq\n'.format(
            d=dataid, g=poly_group, n=poly_name, e=E))
        sys.exit(3)
    sample_subseq = sample_aln_seq[startidx:endidx+1]
    dashcount = sum([x=='-' for x in list(sample_subseq)])
    result = int(dashcount == len(sample_subseq))
    if _debug:
        if result == 1:
            sys.stdout.write('Debug: Sample_DElat_Result {d}, {g}, {n}, E={e}, sample_subseq={v} => {r}\n'.format(
                d=dataid, g=poly_group, n=poly_name, e=E, v=sample_subseq, r=result))
    return result

def Sample_DELany_Result(dataid, poly_group, poly_name, E, I, J):
    'Return 1 if the aligned sample sequence has any deletions in the gene specified.'

    if E['level'] == 'aa':
        if not (dataid, E['gene']) in I or I[(dataid, E['gene'])] is None:
            return 'NA'
        ref_aln_seq = I[(dataid, E['gene'])]['ref']
        sample_aln_seq = I[(dataid, E['gene'])]['sample']
    else:
        if not dataid in J or J[dataid] is None:
            return 'NA'
        ref_aln_seq = J[dataid]['ref']
        sample_aln_seq = J[dataid]['sample']

    result = int('-' in list(sample_aln_seq))
    return result

def Sample_Satisfies_Polymorphism_Criterion(dataid, poly_group, poly_name, E, I, J):
    'Return 1 if the sample satisfies this polymorphism criterion in E, 0 otherwise.'

    result = 0
    if E['type'] == 'SNP':
        result = Sample_SNP_Result(dataid, poly_group, poly_name, E, I, J)
    elif E['type'] == 'DELat':
        result = Sample_DELat_Result(dataid, poly_group, poly_name, E, I, J)
    elif E['type'] == 'DELany':
        result = Sample_DELany_Result(dataid, poly_group, poly_name, E, I, J)
    else:
        sys.stderr.write('Error: Invalid poly_type *{t}*\n'.format(t=E['type']))
        sys.exit(2)
    if _debug:
        sys.stdout.write('Debug: Sample_Satisfies_Polymorphism_Criterion: {d}, {g}, {n} => {r}\n'.format(
            d=dataid, g=poly_group, n=poly_name, r=result))
    return result

def Sample_Has_Polymorphism(dataid, poly_group, poly_name, poly_spec, I, J):
    '''
    If any data for the polymorphism critera were missing, return NA.
    Otherwise, return 1 if all criteria satisfied, otherwise 0.
    '''

    num_criteria = len(poly_spec)
    resultL = []
    for i in range(0, num_criteria):
        resultL.append(Sample_Satisfies_Polymorphism_Criterion(dataid, poly_group, poly_name, poly_spec[i], I, J))

    if 'NA' in resultL:
        result = 'NA'
    else:
        num_criteria_satisfied = sum([v == 1 for v in resultL])
        result = int(num_criteria_satisfied == num_criteria)
    return result

def Compute_SamplePolymorphisms(P, polygroupL, polykeyL, S, dataidL, I, J, aaalndir, ntalndir):
    'Print table of polymorphisms, one row per sample, with boolean values 1/0 (or NA for Not Available).'
    T = []
  # Assemble header row
    header = ['dataid', 'samplename', 'studynumber', 'diagnosis', 'HCC'] + ['('+key[0]+', '+key[1]+')' for key in polykeyL]
    T.append(header)
 # Add one row for each sample
    for dataid in dataidL:
        rowL = [dataid, S[dataid]['samplename'], S[dataid]['studynumber'], S[dataid]['diagnosis'], S[dataid]['hcc']]
        for poly_group, poly_name in polykeyL:
            poly_result = Sample_Has_Polymorphism(dataid, poly_group, poly_name, P[poly_group][poly_name]['spec'], I, J)
            rowL.append(str(poly_result))
        T.append(rowL)
 # Save table to file
    out_path = os.path.join(_outdata_dir, _outdata_prefix+'_sample_polymorphisms.tsv')
    with open(out_path, 'w') as out_fp:
        out_fp.write('\n'.join(['\t'.join([str(x) for x in row]) for row in T])+'\n')
    return T

# Output tables

def SummaryCounts(poly_group, poly_name, T):
    '''
    Given a (group, name) key, count how many samples with HCC=1 vs HCC=0.
    Return { 'HCC': { 'num': n, 'tot': n, 'pct': n }, 'nonHCC': { 'num': n, 'tot': n, 'pct': n } }
    '''
    D = { 'HCC': { 'num': 0, 'tot': 0, 'pct': 0.0 }, 'nonHCC': { 'num': 0, 'tot': 0, 'pct': 0.0 } }

    key = '('+poly_group+', '+poly_name+')'
    j = 0
    while j < len(T[0]) and T[0][j] != key:
        j += 1

    for i in range(1, len(T)):
        rowL = T[i]
        HCC_status = rowL[4]

        if HCC_status == 1:
            if rowL[j] == '1':
                D['HCC']['num'] += 1
            if rowL[j] != 'NA':
                D['HCC']['tot'] += 1

        elif HCC_status == 0:
            if rowL[j] == '1':
                D['nonHCC']['num'] += 1
            if rowL[j] != 'NA':
                D['nonHCC']['tot'] += 1

    D['HCC']['pct'] = round(float(D['HCC']['num']) / D['HCC']['tot']*100.0, _num_dec_pl) if D['HCC']['tot'] else 0.0
    D['nonHCC']['pct'] = round(float(D['nonHCC']['num']) / D['nonHCC']['tot']*100.0, _num_dec_pl) if D['nonHCC']['tot'] else 0.0
    return D

def Compute_SummaryPolymorphisms(T, polykeyL):
    '''
    Return a table of polymorphisms vs HCC and non-HCC counts.
    T = [ [group, name, [HCC_num, HCC_tot, HCC_pct], [nonHCC_num, nonHCC_tot, nonHCC_pct]], ... ]
    '''
    M1 = []
  # Header
    header = ['region', 'polymorphism', 'Freq HCC', 'Freq non-HCC']
    M1.append(header)
  # One row per polymorphism
    for poly_group, poly_name in polykeyL:
        sc = SummaryCounts(poly_group, poly_name, T)
        rowL = [poly_group, poly_name,
            [sc['HCC']['num'], sc['HCC']['tot'], sc['HCC']['pct']],
            [sc['nonHCC']['num'], sc['nonHCC']['tot'], sc['nonHCC']['pct']]]
        M1.append(rowL)
  # Save table to file
    out_path = os.path.join(_outdata_dir, _outdata_prefix+'_summary_polymorphisms.tsv')
    with open(out_path, 'w') as out_fp:
        linecnt = 0
        for rowL in M1:
            linecnt += 1
            if linecnt == 1:
                L = rowL
            else:
                L = [rowL[0], rowL[1],
                    '{a}/{b} ({c:1}%)'.format(a=rowL[2][0], b=rowL[2][1], c=rowL[2][2]),
                    '{a}/{b} ({c:1}%)'.format(a=rowL[3][0], b=rowL[3][1], c=rowL[3][2])]
            out_fp.write('\t'.join([str(x) for x in L])+'\n')
    return M1

def Fishers_Exact_Test(polyyes_nonHCC, polyyes_HCC, polyno_nonHCC, polyno_HCC):
    '''
    H0: polyyes and polyno have the same proportion of nonHCC and HCC samples.
    If p-value < 0.05, reject H0 to conclude they are different.
    Use the two-sided Fisher's exact test, with\
    input table [[polyyes_nonHCC, polyno_nonHCC], [polyyes_HCC, polyno_HCC]]
    '''

    df = pd.DataFrame(
        {'poly-yes':[polyyes_nonHCC, polyyes_HCC], 'poly-no':[polyno_nonHCC, polyno_HCC]},
        index=pd.Index(['no HCC', 'HCC']))

    oddsr, p = fisher_exact(table=df.to_numpy(), alternative='two-sided')
    is_signif = p < _Fishers_Exact_Test_p_value_signif_threshold

    return p, is_signif


def Compute_SummaryPolymorphismsWithPvalues(M1):
    'Append the Fisher\'s exact test p-value for the comparison of the HCC vs non-HCC counts.'
    M2 = []
    header = M1[0] + ['p-value']
    M2.append(header)
    for rowL in M1[1:]:
        region, polymorphism, freq_HCC, freq_nonHCC = rowL
        a, b = freq_HCC[0:2]
        c, d = freq_nonHCC[0:2]
        polyyes_nonHCC = c
        polyyes_HCC = a
        polyno_nonHCC = (d - c) 
        polyno_HCC = (b - a)
        p_value, is_signif = Fishers_Exact_Test(polyyes_nonHCC, polyyes_HCC, polyno_nonHCC, polyno_HCC)
        M2.append(rowL + [round(p_value, 2)])

  # Save table to file
    out_path = os.path.join(_outdata_dir, _outdata_prefix+'_summary_polymorphisms_pvalue.tsv')
    with open(out_path, 'w') as out_fp:
        linecnt = 0
        for rowL in M2:
            linecnt += 1
            if linecnt == 1:
                L = rowL
            else:
                L = [rowL[0], rowL[1],
                    '{a}/{b} ({c:1}%)'.format(a=rowL[2][0], b=rowL[2][1], c=rowL[2][2]),
                    '{a}/{b} ({c:1}%)'.format(a=rowL[3][0], b=rowL[3][1], c=rowL[3][2]),
                    round(rowL[4], 2)]
            out_fp.write('\t'.join([str(x) for x in L])+'\n')
        #for rowL in M2:
        #    out_fp.write('\t'.join([str(x) for x in rowL])+'\n')
    return M2

# Main

if __name__ == '__main__':

  # Parse command-line arguments
    if len(sys.argv) != 10:
        print(_usage_msg)
        sys.exit(1)
    _refname, _refprot_fastafile, _refnucl_fastafile, \
        _polytable_tsvfile, _sampletable_tsvfile, \
        _aaaln_dir, _ntaln_dir, \
        _outdata_dir, _outdata_prefix = sys.argv[1:10]

  # Parse input files
    A, protnameL = Parse_Fasta(_refprot_fastafile, _refname+'_')
    N, dummyN = Parse_Fasta(_refnucl_fastafile)
    P, polygroupL, polykeyL = Parse_PolyTable(_polytable_tsvfile)
    S, sampleL = Parse_SampleTable(_sampletable_tsvfile)
    I = Parse_ProtAlignmentFiles(S, sampleL, protnameL, _aaaln_dir, _refname)
    J = Parse_NtAlignmentFiles(S, _ntaln_dir, _refname)

  # Check that the polymorphism table spec is sensible
    Validate_Polytable_vs_ReferenceSequences(A, P, polygroupL, polykeyL, N, _refname)

  # Tabulate the polymorphisms present in each sample
    T = Compute_SamplePolymorphisms(P, polygroupL, polykeyL, S, sampleL, I, J, _aaaln_dir, _ntaln_dir)

  # Tabulate counts of the polymorphisms present in HCC vs non-HCC samples
    M1 = Compute_SummaryPolymorphisms(T, polykeyL)
    M2 = Compute_SummaryPolymorphismsWithPvalues(M1)


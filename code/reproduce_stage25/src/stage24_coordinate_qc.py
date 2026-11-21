#!/usr/bin/env python3
"""Parse source insert labels; do not infer validated genomic boundaries from sequence offsets."""
import csv,re,json,argparse,zipfile,io
from pathlib import Path
from collections import Counter
ap=argparse.ArgumentParser(description='Audit original GSE110758 source construct labels without projecting allelic offsets')
ap.add_argument('--barcode-map',type=Path,help='Optional uncompressed selected barcode-map CSV')
ap.add_argument('--stage21c-zip',type=Path,default=Path('data/source_archives/HAR_Stage21C_MPRA_HS111_N2_Results.zip'))
ap.add_argument('--outdir',type=Path,default=Path('results/stage24'))
args=ap.parse_args();root=args.outdir;root.mkdir(parents=True,exist_ok=True)
if args.barcode_map:
    stream=io.StringIO(args.barcode_map.read_text(encoding='utf8'))
else:
    with zipfile.ZipFile(args.stage21c_zip) as z:
        stream=io.StringIO(z.read('S1B_Assayed_Barcode_Map.csv').decode('utf8'))
pat=re.compile(r'^(.*?)__species=(human|chimp)__coordinates_tested_in_hg19=\[(chr[^:]+):(\d+)-(\d+)\]__cut=(\d+)_of_(\d+)__length=(\d+)_of_(\d+)__(\d+)-(\d+)__')
insert_counts=Counter()
for x in csv.DictReader(stream):insert_counts[x['insert']]+=1
rows=[]
for ins,nbar in sorted(insert_counts.items()):
 m=pat.match(ins)
 if not m:raise ValueError('Unparsed construct annotation '+ins[:200])
 har,allele,chrom,start,end,cut,ncuts,seqcut,seqtotal,ofs,ofe=m.groups()
 start,end,seqcut,seqtotal,ofs,ofe=map(int,[start,end,seqcut,seqtotal,ofs,ofe])
 if ofe-ofs!=seqcut:raise ValueError('Sequence cut offsets contradict reported length')
 gap=(end-start)-seqtotal
 rows.append(dict(HAR=har,allele=allele,cut=f'{cut}_of_{ncuts}',chrom=chrom,
     reported_hg19_region_start=start,reported_hg19_region_end=end,
     reported_hg19_span=end-start,reported_allele_full_sequence_length=seqtotal,
     reporter_cut_offset_start=ofs,reporter_cut_offset_end=ofe,
     reporter_cut_sequence_length=seqcut,source_design_n_barcodes=nbar,
     hg19_region_span_minus_allele_fullseq_length=gap,
     sequence_offsets_equal_genomic_offsets_if_no_indels=(gap==0),
     possible_indel_or_coordinate_disagreement=(gap!=0),
     projection_status=('SAME_SPAN_ONLY_NOT_GENOMIC_IDENTITY_VALIDATION' if gap==0 else 'SEQUENCE_LENGTH_DIFFERS_DO_NOT_PROJECT_OFFSETS'),
     raw_GEO_construct_annotation=ins))
fields=list(rows[0]);p=root/'S2_255_Source_Construct_AnnotationQC.csv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
summary={'n_inserts':len(rows),'n_HAR':len({x['HAR'] for x in rows}),
         'n_alleles':dict(Counter(x['allele'] for x in rows)),
         'n_span_equal':sum(x['sequence_offsets_equal_genomic_offsets_if_no_indels'] for x in rows),
         'n_span_mismatch':sum(x['possible_indel_or_coordinate_disagreement'] for x in rows),
         'mismatch_HARs':sorted({x['HAR'] for x in rows if x['possible_indel_or_coordinate_disagreement']}),
         'n_mismatch_HARs':len({x['HAR'] for x in rows if x['possible_indel_or_coordinate_disagreement']}),
         'gap_counts':dict(Counter(str(x['hg19_region_span_minus_allele_fullseq_length']) for x in rows)),
         'note':'hg19 source region coordinates and reporter insert allele sequence lengths are distinct measurements; mismatch does not by itself prove a specific indel and sequence offsets are not reliable genomic positions'}
(root/'S2_Annotation_QC_Status.json').write_text(json.dumps(summary,indent=2))
assert len(rows)==255 and summary['n_span_mismatch']==18
print(json.dumps(summary,indent=2))

# 67_scopus_screening_decisions.py
# Record of the title/abstract screening decisions for the Scopus focused search (12 Nov 2026).
# Inputs were the Scopus result list (272 records) and the PubMed screening file (S17b).
# Output: scopus_records_screened.csv (= Supplementary Table S17c). Kept for documentation; paths are those used when it was run.
import json,csv,re
from collections import Counter
recs={int(k):v for k,v in json.load(open('scopus_all.json')).items()}
dup=set(json.load(open('scopus_dups.json')))
pm=list(csv.DictReader(open('/mnt/user-data/outputs/pubmed_focused_records_screened.csv',encoding='utf-8-sig')))
norm=lambda s: re.sub(r'[^a-z0-9]','',s.lower())[:60]
pmd={norm(r['title']):r for r in pm}
LB ="Exclude (E3): liquid biopsy / circulating tumour cells / extracellular vesicles / circulating nucleic acids - not a tumour-microenvironment chip"
MEAS="Exclude (I2): microfluidics used only as a measurement technology on patient samples or cell lines - not an organ-on-chip model"
DEV="Exclude (E2): device or single-cell method (isolation, culture, sequencing, mechanotyping, sensing) without a tumour-microenvironment readout"
NOCHIP="Exclude (I2): not a microfluidic/organ-on-chip model study"
BG="Background (I1): review, editorial, book chapter or survey - not primary research"
COR="Exclude (I1): correction notice"
OTH="Exclude (E1): not colorectal cancer"
ELIG_LOW="Eligible (low priority): CRC chip/microfluidic organoid study without tumour-microenvironment co-culture"
ELIG="Eligible: CRC chip or microphysiological study"
D={}
def s(lst,dec,why=''):
    for i in lst: D[i]=(dec,why)
s([3,31,32,50,75,80,83,98,99,149,152,155,162,191,206,207,235,240,254,257,258,261,265,268],LB)
s([13,14,26,46,62,70,82,104,112,114,203,238,245,249,259,269],MEAS)
s([12,72,93,116,154,189,210,216,220,223,224,230,233,248,253],DEV)
s([8,9,10,15,22,27,29,35,49,52,53,61,68,74,76,77,79,86,88,96,97,100,101,102,106,108,109,110,111,113,119,120,121,123,125,128,130,132,133,134,135,137,138,139,140,141,142,143,145,146,147,150,151,153,157,159,160,163,166,167,168,171,173,175,176,178,181,183,184,186,187,188,192,194,196,199,201,208,212,218,219,221,225,227,228,241,246,251,252,263,266,270],BG)
s([6,73],COR)
s([4],BG)
s([92],ELIG,'Colonocyte-microbiota co-culture chip; no patient-derived marker selection')
s([54,55,78,87,156],OTH)
s([16],NOCHIP,'Immunotherapy mechanism study (abstract screen)')
s([19],NOCHIP,'Stool-collection feasibility study')
s([24],NOCHIP,'Organoid drug-sensitivity profiling without a chip; liver metastasis')
s([71],NOCHIP,'Macrophage-derived nanovesicle therapy study')
s([84],NOCHIP,'Drug tolerance across bacteria, fungi and cancer cells')
s([91],NOCHIP,'Murine single-cell RNA-seq; no chip')
s([115],NOCHIP,'Computational ligand-receptor method for bulk and spatial transcriptomics')
s([126],NOCHIP,'Patient monocyte programming study')
s([136],NOCHIP,'Patient CRC single-cell atlas (CMS refinement); no chip - relevant to patient-data background only')
s([169],NOCHIP,'MSI/MMR testing in solid tumours')
s([11,58],ELIG_LOW)
s([2],ELIG,'Colon tumour spheroids in magnetic microcontainers on chip; no patient-derived marker selection')
s([57],ELIG,'GAD1 pathway expression linked to survival in trial patients (CALGB/SWOG 80405) and studied in a CRC tumour-chip invasion model. Already in S17')
s([63],ELIG,'Protocol for immunocompetent tumour-on-chip from patient tumours (autologous tumour cells and CD8+ TILs); no reference marker panel')
s([95],ELIG,'Patient-derived tumour spheroid-induced angiogenesis platform; no patient-derived marker panel')
s([172],ELIG,'Vascularised micro-tumour (CRC) recapitulating in vivo drug responses; validation against xenografts, not patient transcriptomic programs. Already in S17')
s([174],ELIG,'Pumpless multi-organ MPS for metabolism-dependent chemotherapy toxicity; no TME marker readout')
only=[n for n in sorted(recs) if n not in dup]
assert sorted(D)==only, (set(only)-set(D), set(D)-set(only))
GAP={57:'Partial: pathway genes from patient trial data examined on chip; no compartment panel, no independent single-cell or spatial validation',
     172:'Partial: chip validated against in vivo drug response, not against patient transcriptomic programs'}
out=[]
for n in sorted(recs):
    r=recs[n]; o=dict(scopus_no=n,year=r['year'],first_author=r['first_author'],title=r['title'],source=r['source'],doc_type=r['type'])
    if n in dup:
        p=pmd[norm(r['title'])]
        o.update(in_pubmed_set='Yes',pmid=p['pmid'],decision='Duplicate of PubMed record',reason='Screened in S17b: '+p['decision'],gap_relevance=p['gap_relevance'])
    else:
        dec,why=D[n]
        o.update(in_pubmed_set='No',pmid='',decision=dec.split(':')[0],reason=dec.split(': ',1)[1]+('; '+why if why else ''),
                 gap_relevance=GAP.get(n,'No: does not derive patient-data-based marker panels for chip assessment'))
    out.append(o)
cols=['scopus_no','year','first_author','title','source','doc_type','in_pubmed_set','pmid','decision','reason','gap_relevance']
with open('/mnt/user-data/outputs/scopus_records_screened.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(out)
new=[o for o in out if o['in_pubmed_set']=='No']
c=Counter(o['decision'] for o in new); print(len(out),len(new)); print(c)
print('elig new', sum(o['decision'].startswith('Eligible') for o in new), 'partial', sum(o['gap_relevance'].startswith('Partial') for o in new))

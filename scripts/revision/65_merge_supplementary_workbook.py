"""65_merge_supplementary_workbook.py - build the revised Supplementary Tables workbook (V5) with openpyxl.

Used instead of an R/openxlsx script because openxlsx corrupts this workbook on re-save
(it writes links to drawing files that it does not include).
  1. corrects TNF/NF-kB 0.38-0.49 -> 0.38-0.50 in every sheet (S11, S12)
  2. adds the revision columns to S11 and S12 from the revised evidence tables
  3. adds the revision sheets from script 64 and the input manifest from script 63
     (random panels = S15, input manifest = S16; S17 reserved for the prior-art table)
  4. sorts sheets numerically and rebuilds Contents with hyperlinks
Run from the project root (Python 3 with openpyxl:  pip install openpyxl):
  python 05_analysis_scripts/revision/65_merge_supplementary_workbook.py ^
    10_submission/Supplementary_Tables_aem.xlsx ^
    08_tables/revision_2026-11/Supplementary_Tables_revision.xlsx ^
    06_results/revision_2026-11/63_input_manifest_sha256.csv ^
    08_tables/Table_program_evidence_8.csv 08_tables/Table_chip_assay_proposal_8.csv ^
    10_submission/Supplementary_Tables_V5_aem.xlsx
"""
import sys, csv, re, copy, openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
BASE, REV, MAN, EVT, CHT, OUT = sys.argv[1:7]
import zipfile, io
def safe_load(path, **kw):
    """Open an .xlsx; if it carries links to drawing files that are missing (an openxlsx bug), drop those links first."""
    try:
        return openpyxl.load_workbook(path, **kw)
    except KeyError:
        zin = zipfile.ZipFile(path); names = set(zin.namelist()); buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zout:
            for it in zin.infolist():
                d = zin.read(it.filename)
                if it.filename.endswith('.rels'):
                    t = d.decode('utf8'); base = it.filename.rsplit('/_rels/', 1)[0]
                    def ok(m):
                        tgt = m.group(1)
                        full = tgt.lstrip('/') if tgt.startswith('/') else (base + '/' + tgt)
                        parts = []
                        for seg in full.split('/'):
                            if seg == '..': parts.pop()
                            elif seg != '.': parts.append(seg)
                        return '/'.join(parts) in names or 'TargetMode="External"' in m.group(0)
                    t = re.sub(r'<Relationship [^>]*Target="([^"]+)"[^>]*/>', lambda m: m.group(0) if ok(m) else '', t)
                    d = t.encode('utf8')
                elif it.filename.startswith('xl/worksheets/sheet') and it.filename.endswith('.xml'):
                    t = d.decode('utf8')
                    rel = it.filename.replace('worksheets/', 'worksheets/_rels/') + '.rels'
                    if rel in names:
                        rt = zin.read(rel).decode('utf8')
                        for rid, tgt in re.findall(r'Id="([^"]+)"[^>]*Target="([^"]+)"', rt):
                            if 'drawing' in tgt and ('xl/drawings/' + tgt.split('/')[-1]) not in names:
                                t = re.sub(r'<(legacyDrawing|drawing) r:id="%s"/>' % re.escape(rid), '', t)
                    d = t.encode('utf8')
                zout.writestr(it.filename, d)
        buf.seek(0); print('NOTE: repaired broken drawing links in', path)
        return openpyxl.load_workbook(buf, **kw)
wb = safe_load(BASE)
HF = Font(bold=True); HFILL = PatternFill('solid', fgColor='D9E2F3'); WRAP = Alignment(wrap_text=True, vertical='top')
BOT = Border(bottom=Side(style='thin'))
log = []
# 1. 0.49 -> 0.50
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str) and ('0.38-0.49' in c.value or '0.38–0.49' in c.value):
                c.value = c.value.replace('0.38-0.49', '0.38-0.50').replace('0.38–0.49', '0.38–0.50'); log.append(f'fixed {ws.title}!{c.coordinate}')
def rows(path): return list(csv.DictReader(open(path, encoding='utf-8-sig')))
def header_row(ws, first):
    for r in range(1, 20):
        if ws.cell(r, 1).value == first: return r
    raise SystemExit(f'header {first} not found in {ws.title}')
def style_like(src, dst):
    dst.font = copy.copy(src.font); dst.fill = copy.copy(src.fill); dst.border = copy.copy(src.border); dst.alignment = copy.copy(src.alignment)
# 2. S11 new columns
ev = {r['program']: r for r in rows(EVT)}
ws = wb[[n for n in wb.sheetnames if n.startswith('S11')][0]]; hr = header_row(ws, 'Program code'); nc = ws.max_column
add = [('External cohorts (k)','k_cohorts',10), ('Validation patients (n)','n_validation_patients',10), ('Random-panel p (script 53)','random_panel_p',11),
       ('Random-panel comparison','random_panel_detail',34), ('Patients positive, criterion A','spatial_A_positive',11), ('Patients positive, criterion B','spatial_B_positive',11),
       ('Spatial label, random-selection reference','spatial_label_random_reference',26), ('Evidence qualifier','evidence_qualifier',60)]
for j, (h, k, w) in enumerate(add, start=nc + 1):
    c = ws.cell(hr, j, h); style_like(ws.cell(hr, nc), c); ws.column_dimensions[get_column_letter(j)].width = w
    for r in range(hr + 1, ws.max_row + 1):
        p = ws.cell(r, 1).value
        if p in ev:
            v = ev[p][k]; v = int(v) if k in ('k_cohorts', 'n_validation_patients') else (float(v) if k == 'random_panel_p' else v)
            c = ws.cell(r, j, v); style_like(ws.cell(r, nc), c); c.alignment = WRAP
log.append(f'S11: {len(add)} columns added')
# 3. S12 updates
ch = {r['gene']: r for r in rows(CHT)}
ws = wb[[n for n in wb.sheetnames if n.startswith('S12')][0]]; hr = header_row(ws, 'Gene'); nc = ws.max_column
hdr = {ws.cell(hr, j).value: j for j in range(1, nc + 1)}
c = ws.cell(hr, nc + 1, 'Evidence qualifier'); style_like(ws.cell(hr, nc), c); ws.column_dimensions[get_column_letter(nc + 1)].width = 40
n = 0
for r in range(hr + 1, ws.max_row + 1):
    g = ws.cell(r, 1).value
    if g in ch:
        ws.cell(r, hdr['Spatial evidence']).value = ch[g]['spatial_evidence']
        ws.cell(r, hdr['Interpretation limits']).value = ch[g]['interpretation_limits']
        c = ws.cell(r, nc + 1, ch[g]['evidence_qualifier']); style_like(ws.cell(r, nc), c); c.alignment = WRAP; n += 1
log.append(f'S12: {n} genes updated')
# 4. revision sheets
TITLE = {
 'S1d_patients_by_program': 'Supplementary Table S1d. Independent validation patients by program and cohort',
 'S7c_alternative_panels': 'Supplementary Table S7c. Alternative panels from frequently selected substitutes (post hoc)',
 'S7d_mmr_site_strata': 'Supplementary Table S7d. Discovery cross-validation by mismatch-repair status and tumour site (post hoc, uncorrected)',
 'S8c_calibration_transport': 'Supplementary Table S8c. Calibration and transport of the frozen scale by cohort (post hoc)',
 'S8d_leave_one_donor_out': 'Supplementary Table S8d. Leave-one-donor-out and jackknife analysis (post hoc)',
 'S8e_meta_sensitivity': 'Supplementary Table S8e. Meta-analysis sensitivity: REML, Hartung–Knapp, Paule–Mandel, leave-one-cohort-out (post hoc)',
 'S8f_multiplicity': 'Supplementary Table S8f. Multiplicity adjustment across program × cohort tests (post hoc)',
 'S10e_spatial_table_sens': 'Supplementary Table S10e. Spatial thresholds and cohort exclusion (post hoc)',
 'S10f_spatial_rerun': 'Supplementary Table S10f. Spatial re-analysis: block size and residualised model (post hoc)',
 'S10g_spatial_negcontrol': 'Supplementary Table S10g. Spatial negative-control targets (post hoc)',
 'S14_amendment_log': 'Supplementary Table S14. Protocol changes made after related results were known',
 'S15_null_random_panels': 'Supplementary Table S15. Comparison with expression-matched random panels and negative-control targets (post hoc)',
 'S16_input_manifest': 'Supplementary Table S16. Input manifest: SHA-256 checksums of original data files'}
tfont = copy.copy(wb[wb.sheetnames[1]]['A1'].font)
def add_sheet(name, header, data):
    if name in wb.sheetnames: del wb[name]
    ws = wb.create_sheet(name); ws['A1'] = TITLE.get(name, name); ws['A1'].font = tfont
    for j, h in enumerate(header, 1):
        c = ws.cell(3, j, h); c.font = HF; c.fill = HFILL; c.border = BOT; c.alignment = WRAP
        ws.column_dimensions[get_column_letter(j)].width = min(max(len(str(h)) + 2, 10), 40)
    for i, row in enumerate(data, 4):
        for j, v in enumerate(row, 1): ws.cell(i, j, v)
    ws.freeze_panes = 'A4'; log.append(f'added {name} ({len(data)} rows)')
rv = safe_load(REV, data_only=True)
for s in rv.worksheets:
    vals = [list(r) for r in s.iter_rows(values_only=True)]
    vals = [r for r in vals if any(v is not None for v in r)]
    if vals: add_sheet({'S16_null_random_panels': 'S15_null_random_panels'}.get(s.title, s.title), vals[0], vals[1:])
if MAN != '-':
    m = list(csv.reader(open(MAN, encoding='utf-8-sig')))
    add_sheet('S16_input_manifest', m[0], m[1:])
# 5. order + Contents
def key(s):
    if s == 'Contents': return -1
    if s == 'References': return 1e6
    mm = re.match(r'^S(\d+)([a-z]?)_', s)
    return int(mm.group(1)) * 100 + (ord(mm.group(2)) - 96 if mm.group(2) else 0) if mm else 9e5
wb._sheets.sort(key=lambda w: key(w.title))
ct = wb['Contents']; hr = header_row(ct, 'Table')
abb_r = next((r for r in range(hr, ct.max_row + 1) if str(ct.cell(r, 1).value or '').startswith('Abbreviations')), None)
abb = [(ct.cell(r, 1).value, copy.copy(ct.cell(r, 1).font)) for r in range(abb_r, ct.max_row + 1)] if abb_r else []
abb = [a for a in abb if a[0]]
sty = [copy.copy(ct.cell(hr + 1, j)._style) for j in range(1, 5)]
link_font = copy.copy(ct.cell(hr + 1, 2).font)
for r in range(hr + 1, ct.max_row + 60):
    for j in range(1, 6): ct.cell(r, j).value = None; ct.cell(r, j).hyperlink = None
r = hr + 1
for s in wb.sheetnames:
    if s == 'Contents': continue
    w = wb[s]; t = str(w['A1'].value or s)
    lab = 'References' if s == 'References' else re.sub(r'^(Supplementary Table S\d+[a-z]?)\..*$', r'\1', t)
    short = 'Works cited in these tables' if s == 'References' else re.sub(r'^Supplementary Table S\d+[a-z]?\.\s*', '', t)
    hrow = next((i for i in range(2, 20) if w.cell(i, 1).value is not None and i > 1 and w.cell(i, 2).value is not None), 3)
    nrows = max(0, sum(1 for rr in w.iter_rows(min_row=hrow + 1, values_only=True) if any(v is not None for v in rr)))
    for j, v in enumerate([lab, f'=HYPERLINK("#\'{s}\'!A1","{s}")', short, nrows], 1):
        c = ct.cell(r, j, v); c._style = copy.copy(sty[j - 1])
    r += 1
r += 1
for v, f in abb:
    c = ct.cell(r, 1, v); c.font = f; r += 1
ct.cell(r + 1, 1, 'Post hoc: analyses added in response to pre-submission review (seed 20261111), specified after the primary results were known.')
wb.active = 0
wb.save(OUT)
print('\n'.join(log)); print('sheets:', len(wb.sheetnames)); print(wb.sheetnames)

"""Reproducible Girard T model reaction and schematic steroid MS³ workflow."""
import json,math
from pathlib import Path
from collections import Counter
from dataclasses import replace
from xml.etree import ElementTree as ET
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from codex_organic_chem.figure_style import render_molecule,resolve_preset,_strip_to_body
from codex_organic_chem.route_figure import _arrow_svg
from codex_organic_chem.figure_audit import FigureElement,audit_figure,svg_to_png,raster_ink_check

BLUE=(0.03,0.18,0.65);RED=(0.75,0.04,0.06);DARK=(0.16,0.16,0.16)
preset=resolve_preset('acs',overrides={'bond_length_px':36,'max_font_size':18,'min_font_size':11})
preset['explicit_methyl']=True
from rdkit.Chem import rdMolDescriptors
from html import escape
ORANGE=(0.8,0.36,0.02)
strings=['[CH3:1][C:2](=[O:3])[CH3:4]','[NH2:5][NH:6][C:7](=[O:8])[CH2:9][N+:10]([CH3:11])([CH3:12])[CH3:13]','[CH3:1][C:2](=[N:5][NH:6][C:7](=[O:8])[CH2:9][N+:10]([CH3:11])([CH3:12])[CH3:13])[CH3:4]','[OH2:3]','[CH3:1][C:2]([OH:3])([NH:5][NH:6][C:7](=[O:8])[CH2:9][N+:10]([CH3:11])([CH3:12])[CH3:13])[CH3:4]']
ms=[Chem.MolFromSmiles(s) for s in strings]
def inventory(mols):
 c=Counter()
 for m in mols:
  for a in m.GetAtoms():c[a.GetSymbol()]+=1;c['H']+=a.GetTotalNumHs()
 return c,sum(Chem.GetFormalCharge(m) for m in mols)
assert inventory(ms[:2])==inventory(ms[2:4])==inventory([ms[4]])
for m in ms: print(rdMolDescriptors.CalcMolFormula(m),Chem.GetFormalCharge(m))
for m in (ms[1],ms[2],ms[4]):
 ids={a.GetAtomMapNum():a for a in m.GetAtoms()}
 assert ids[10].GetFormalCharge()==1 and ids[10].GetDegree()==4
 assert m.GetBondBetweenAtoms(ids[7].GetIdx(),ids[8].GetIdx()).GetBondTypeAsDouble()==2
assert ms[4].GetAtomWithIdx(1).GetDegree()==4
preset['max_font_size']=16

coords=[{2:(0,0),3:(0,1),1:(-.866,-.5),4:(.866,-.5)}, {11:(-1,0),10:(0,0),12:(0,1),13:(0,-1),9:(1,0),7:(2,0),8:(2,1),6:(3,0),5:(4,0)}, {1:(-1,0),2:(0,0),4:(0,1),5:(1,0),6:(2,0),7:(3,0),8:(3,1),9:(4,0),10:(5,0),11:(6,0),12:(5,1),13:(5,-1)}, {3:(0,0)}]
def draw(m,xy):
 m=Chem.Mol(m);ids={a.GetAtomMapNum():a.GetIdx() for a in m.GetAtoms()}
 conf=Chem.Conformer(m.GetNumAtoms());conf.Set3D(False)
 for n,i in ids.items():
  conf.SetAtomPosition(i,(*xy[n],0));a=m.GetAtomWithIdx(i)
  if n in (2,7):a.SetProp('atomLabel','C')
  if n==9:a.SetProp('atomLabel','CH<sub>2</sub>')
  if n==6:a.SetProp('atomLabel','NH')
  if n==5 and a.GetTotalNumHs()==2:a.SetProp('atomLabel','NH<sub>2</sub>')
  a.SetAtomMapNum(0)
 m.RemoveAllConformers();m.AddConformer(conf)
 return render_molecule(mol=m,preset=preset),ids
items=[draw(m,c) for m,c in zip(ms[:4],coords)]
a,g,h,w=[d for d,ids in items]

# Run from repository root: .venv/bin/python figures/girard_t_isomer_A.py
from codex_organic_chem.service import chem_normalize_structure,chem_reaction_analyze
for smi in strings[:4]:
    assert chem_normalize_structure(smi).get('canonical_smiles')
assert chem_reaction_analyze('.'.join(strings[:2])+'>>'+'.'.join(strings[2:4]))['analysis']['mass_charge_balance']['status']=='balanced'
width=1100;height=494;y=108
parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
'<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}.head{font-size:16px;font-weight:600;fill:#294157}.stage{font-size:15px;font-weight:600;fill:#294157}.resultA{font-size:17px;font-weight:600;fill:#3879B0}.resultB{font-size:17px;font-weight:600;fill:#CC823A}.label{font-size:14px}.small{font-size:12px}.tiny{font-size:11px}.blue{fill:#234B79}.rxnArrow,.arrowHead{stroke:#111;stroke-width:1.6;fill:none}.arrowLabel{font-size:13px}.axis{stroke:#444;stroke-width:1;fill:none}.peak{stroke-width:2;fill:none}</style>',f'<rect width="{width}" height="{height}" fill="white"/>']
els=[]
def text(x,y,s,cls='label',anchor='middle'):
    parts.append(f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{escape(s)}</text>')
def place(i,x,anchor,label):
    d,ids=items[i];yy=y-d.atom_points[ids[anchor]][1]
    parts.append(f'<g transform="translate({x},{yy})">{d.svg_body}</g>')
    geo=d.translated_geometry(x,yy)
    els.append(FigureElement(kind='molecule',label=label,box=geo['ink_box'],bond_lengths_px=geo['bond_lengths_px'],label_boxes=geo['label_boxes']))
def arrow(x1,yy,x2):
    parts.append(f'<path d="M {x1} {yy} H {x2-7}" stroke="#71818B" stroke-width="1.5" fill="none"/><path d="M {x2-7} {yy-3.5} L {x2} {yy} L {x2-7} {yy+3.5} Z" fill="#71818B"/>')
def box(x,yy,w,h,label,kind='input'):
    color={'input':'#C5D2DD','A':'#3879B0','B':'#CC823A'}[kind]
    fill='#F5F7F9' if kind=='input' else '#FFFFFF'
    parts.append(f'<rect x="{x}" y="{yy}" width="{w}" height="{h}" fill="{fill}" stroke="{color}" stroke-width="1.2"/>')
    text(x+w/2,yy+h/2+5,label,'stage' if kind=='input' else 'result'+kind)
def plus(xx):text(xx,y+6,'+')
text(26,30,'A','head','start');text(width/2,30,'Girard T derivatization (model reaction)','head')
rowwidth=a.width+g.width+h.width+w.width+2*44+108
x=(width-rowwidth)/2
place(0,x,2,'acetone model');text(x+a.width/2,181,'Acetone','small');x+=a.width+22;plus(x);x+=22
place(1,x,10,'Girard T cation');text(x+g.width/2,181,'Girard T reagent','small');x+=g.width+18
arrow(x,y,x+72);x+=90
hx=x;place(2,x,2,'Girard T hydrazone');text(x+h.width/2,181,'Charge-tagged hydrazone','small')
idx=items[2][1][10];nx=hx+h.atom_points[idx][0]
text(nx,54,'Permanent positive charge','small')
x+=h.width+22;plus(x);x+=22;place(3,x,3,'water')
text(width/2,204,'Acetone: model carbonyl compound; analytical workflow: steroid isomers.','small')
# Restrained thin-box styling following the user's reference.
text(170,247,'Enhanced ionization','stage')
text(565,247,'MS³ fragmentation','stage')
text(932,247,'Isomer-specific identification','stage')
arrow(351,350,401)
def spectrum(x,baseline,w,h,peaks,color,axes=True,xlabel=True):
    if axes:parts.append(f'<path class="axis" d="M {x} {baseline-h} V {baseline} H {x+w}"/>')
    else:parts.append(f'<path class="axis" d="M {x} {baseline} H {x+w}"/>')
    for pos,amp in peaks:
        parts.append(f'<path class="peak" stroke="{color}" d="M {x+pos*w:.2f} {baseline} V {baseline-amp*h:.2f}"/>')
    if xlabel:text(x+w/2,baseline+18,'m/z','tiny')
# Non-experimental example targets shift in m/z upon derivatization.
pre=[(.16,.08),(.31,.23),(.60,.10),(.84,.07)]
post=[(.16,.08),(.41,.10),(.69,.85),(.89,.08)]
spectrum(40,382,120,84,pre,'#8995A0');spectrum(199,382,120,84,post,'#3879B0')
text(100,284,'Before derivatization','small');text(259,284,'After derivatization','small');arrow(169,343,190)
text(170,423,'Charge tagging shifts the target mass','tiny')
text(170,447,'Schematic spectra (not experimental data)','tiny')
box(441,270,266,35,'Derivatized steroid ion')
text(510,332,'MS²','small');arrow(535,328,576);text(602,332,'MS³','small')
patternA=[(.12,.18),(.30,.75),(.56,.29),(.79,.49)]
patternB=[(.12,.18),(.41,.48),(.64,.78),(.89,.29)]
text(436,372,'A','small');spectrum(461,387,265,42,patternA,'#3879B0',False,False)
text(436,424,'B','small');spectrum(461,439,265,42,patternB,'#CC823A',False,False)
text(572,459,'Schematic MS³ fragmentation patterns','tiny')
# Diagnostic pattern to identity: two separate arrows, no mass-based separation.
arrow(745,368,812);arrow(745,420,812)
box(835,346,200,43,'Isomer A','A');box(835,398,200,43,'Isomer B','B')
text(935,296,'Same formula','small');text(935,314,'Same precursor m/z','small')
text(935,459,'Structure-dependent fragmentation','tiny')
text(25,477,'Conceptual workflow • Takeo et al., Anal. Chem. 2019 • DOI: 10.1021/acs.analchem.9b00619','tiny','start')
parts.append('</svg>')
root=ET.fromstring('\n'.join(parts));meta=ET.SubElement(root,'{http://www.w3.org/2000/svg}metadata')
meta.text=json.dumps({'mapped_reaction':strings[:4],'element_hydrogen_charge_balance':'balanced; +1 on both sides','model':'acetone only; not actual steroid analyte','isomers':'Conceptual A and B are defined to share molecular formula, derivatization stoichiometry and precursor charge; no particular steroid or fragment assignment','spectra':'Entirely schematic; normalized arbitrary display coordinates, not experimental intensities or m/z; no enhancement factor asserted','reference':{'doi':'10.1021/acs.analchem.9b00619','experimental_fragment_structures':'not obtained; none drawn'},'source_file':'girard_t_isomer_A.py'},ensure_ascii=False)
ET.register_namespace('','http://www.w3.org/2000/svg')
out=Path(__file__).resolve().with_suffix('.svg');out.write_text(ET.tostring(root,encoding='unicode'));png=out.with_suffix('.png')
assert not svg_to_png(out,png,scale=4800/width)
with Image.open(png) as im:im.load();im.save(png,dpi=(600,600))
with Image.open(png) as im:im.load();assert im.format=='PNG';print('PNG',im.size,im.info['dpi'])
warn,ink=raster_ink_check(png);assert not warn and ink['status']=='ok'
report=audit_figure(renderer='rdkit_Girard_T_A',width=width,height=height,elements=els,target_bond_px=36,png_path=png)
assert not report['blocking_issues'] and not report['warnings'],report
assert all(p.stat().st_size>0 for p in (out,png));ET.parse(out)
print('Verified:',out,ink)

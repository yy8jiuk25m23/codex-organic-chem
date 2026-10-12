"""Reproduce: .venv/bin/python figures/fdaa_tims_alanine.py (from repository root).
Chemistry: user-supplied stereospecific SMILES, independently checked by RDKit.
Literature figures unavailable through current web proxy. Mobility trace schematic.
"""
import json,math
from pathlib import Path
from dataclasses import replace
from collections import Counter
from html import escape
from xml.etree import ElementTree as ET
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors,rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from codex_organic_chem.service import chem_normalize_structure,chem_reaction_analyze
from codex_organic_chem.figure_style import prepare_mol,render_molecule,resolve_preset,_strip_to_body
from codex_organic_chem.figure_audit import FigureElement,audit_figure,svg_to_png,raster_ink_check
S={
'L-FDAA':'C[C@@H](C(=O)N)NC1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])F',
'D-Ala':'N[C@H](C)C(=O)O',
'L-Ala':'N[C@@H](C)C(=O)O',
'L-FDAA-D-Ala':'C[C@@H](C(=O)N)NC1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])N[C@H](C)C(=O)O',
'L-FDAA-L-Ala':'C[C@@H](C(=O)N)NC1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])N[C@@H](C)C(=O)O'}
M={k:Chem.MolFromSmiles(v) for k,v in S.items()}
expected={'L-FDAA':['S'],'D-Ala':['R'],'L-Ala':['S'],'L-FDAA-D-Ala':['S','R'],'L-FDAA-L-Ala':['S','S']}
for k,m in M.items():
 assert chem_normalize_structure(S[k]).get('canonical_smiles')
 Chem.AssignStereochemistry(m,cleanIt=True,force=True)
 assert [c for i,c in Chem.FindMolChiralCenters(m)]==expected[k]
 assert Chem.GetFormalCharge(m)==0
 print(k,rdMolDescriptors.CalcMolFormula(m),Chem.FindMolChiralCenters(m))
f=M['L-FDAA'];fi=next(a.GetIdx() for a in f.GetAtoms() if a.GetSymbol()=='F');ar=f.GetAtomWithIdx(fi).GetNeighbors()[0].GetIdx()
assert ar==8 and f.GetAtomWithIdx(ar).GetIsAromatic()
# Original ring positions: NH-alaninamide at 1, nitro at 2/4, F at 5.
ring=f.GetRingInfo().AtomRings()[0];assert len(ring)==6
for aa,p in [('D-Ala','L-FDAA-D-Ala'),('L-Ala','L-FDAA-L-Ala')]:
 m=M[p];assert not any(a.GetSymbol()=='F' for a in m.GetAtoms())
 assert rdMolDescriptors.CalcMolFormula(m)=='C12H15N5O7'
 assert m.GetBondBetweenAtoms(ar,18).GetBondTypeAsDouble()==1
 assert m.GetAtomWithIdx(18).GetSymbol()=='N' and m.GetAtomWithIdx(18).GetTotalNumHs()==1
 assert m.GetBondBetweenAtoms(18,19) and m.GetAtomWithIdx(19).HasProp('_CIPCode')
 for i in range(18):
  assert f.GetAtomWithIdx(i).GetAtomicNum()==m.GetAtomWithIdx(i).GetAtomicNum()
 for b in f.GetBonds():
  i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
  if fi not in (i,j):assert m.GetBondBetweenAtoms(i,j).GetBondType()==b.GetBondType()
 assert chem_reaction_analyze(S['L-FDAA']+'.'+S[aa]+'>>'+S[p]+'.F')['analysis']['mass_charge_balance']['status']=='balanced'
assert rdMolDescriptors.CalcExactMolWt(M['L-FDAA-D-Ala'])==rdMolDescriptors.CalcExactMolWt(M['L-FDAA-L-Ala'])
P=resolve_preset('acs',overrides={'bond_length_px':30,'min_font_size':10,'max_font_size':14,'bond_line_width':1.4});P['explicit_methyl']=True
# Shared fixed coordinates keep both reaction rows aligned; RDKit computes each wedge.
for m in M.values():rdDepictor.Compute2DCoords(m)
ref=M['L-FDAA-D-Ala']
rdDepictor.GenerateDepictionMatching2DStructure(ref,f,atomMap=[(i,i) for i in range(18)])
M['L-FDAA-L-Ala'].RemoveAllConformers();M['L-FDAA-L-Ala'].AddConformer(Chem.Conformer(ref.GetConformer()))
M['L-Ala'].RemoveAllConformers();M['L-Ala'].AddConformer(Chem.Conformer(M['D-Ala'].GetConformer()))
# Orient the product's principal axis horizontally, applying the same rotation to FDAA.
xy=[(ref.GetConformer().GetAtomPosition(i).x,ref.GetConformer().GetAtomPosition(i).y) for i in range(ref.GetNumAtoms())]
mx=sum(x for x,y in xy)/len(xy);my=sum(y for x,y in xy)/len(xy)
xx=sum((x-mx)**2 for x,y in xy);yy=sum((y-my)**2 for x,y in xy);cross=sum((x-mx)*(y-my) for x,y in xy)
angle=-.5*math.atan2(2*cross,xx-yy);co=math.cos(angle);si=math.sin(angle)
for key in ('L-FDAA','L-FDAA-D-Ala','L-FDAA-L-Ala'):
 conf=M[key].GetConformer()
 for i in range(M[key].GetNumAtoms()):
  v=conf.GetAtomPosition(i);conf.SetAtomPosition(i,(co*v.x-si*v.y,si*v.x+co*v.y,0))
BLUE=(.18,.43,.67);ORANGE=(.8,.45,.15)
def draw(k):
 m=Chem.Mol(M[k]);conf=m.GetConformer();ls=sorted((conf.GetAtomPosition(b.GetBeginAtomIdx())-conf.GetAtomPosition(b.GetEndAtomIdx())).Length() for b in m.GetBonds());unit=ls[len(ls)//2]
 for i in range(m.GetNumAtoms()):
  v=conf.GetAtomPosition(i);conf.SetAtomPosition(i,(v.x/unit,v.y/unit,0))
 idx=1 if k in ('D-Ala','L-Ala') else 19 if k.startswith('L-FDAA-') else None
 if idx is not None:m.GetAtomWithIdx(idx).SetProp('atomLabel','CH')
 d=render_molecule(mol=m,preset=P)
 if idx is None:return d
 colors={idx:BLUE if 'D-Ala' in k else ORANGE};bc={b.GetIdx():(0,0,0) for b in m.GetBonds()}
 dr=rdMolDraw2D.MolDraw2DSVG(d.width,d.height);o=dr.drawOptions()
 for n,v in {'fixedBondLength':30,'bondLineWidth':1.4,'baseFontSize':P['base_font_size'],'minFontSize':10,'maxFontSize':14,'padding':0,'centreMoleculesBeforeDrawing':True,'clearBackground':False,'scaleBondWidth':False,'explicitMethyl':True,'circleAtoms':False,'continuousHighlight':False,'fillHighlights':False,'highlightBondWidthMultiplier':8,'scaleHighlightBondWidth':False}.items():setattr(o,n,v)
 o.useBWAtomPalette();rdMolDraw2D.PrepareAndDrawMolecule(dr,m,highlightAtoms=list(colors),highlightBonds=list(bc),highlightAtomColors=colors,highlightBondColors=bc);dr.FinishDrawing()
 return replace(d,svg_body=_strip_to_body(dr.GetDrawingText()))
D={k:draw(k) for k in M}
aw=max(D['D-Ala'].width,D['L-Ala'].width);fw=D['L-FDAA'].width;pw=max(D['L-FDAA-D-Ala'].width,D['L-FDAA-L-Ala'].width)
fw+=95;pw+=95
rowh=max(d.height for d in D.values())+38
width=max(1100,math.ceil(aw+fw+pw+240));height=math.ceil(92+2*rowh+216)
parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<style>text{font-family:Arial,Helvetica,sans-serif;fill:#151515}.title{font-size:17px;font-weight:600;fill:#294157}.heading{font-size:15px;font-weight:600;fill:#294157}.label{font-size:13px}.small{font-size:11px}.blue{fill:#2E6EAB}.orange{fill:#CC7326}.axis{stroke:#555;stroke-width:1.1;fill:none}</style>',f'<rect width="{width}" height="{height}" fill="white"/>'];elements=[]
def text(x,y,s,cls='label',anchor='middle'):parts.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" class="{cls}">{escape(s)}</text>')
def arrow(x1,y,x2):parts.append(f'<path d="M {x1} {y} H {x2-8}" stroke="#71818B" stroke-width="1.5" fill="none"/><path d="M {x2-8} {y-4} L {x2} {y} L {x2-8} {y+4} Z" fill="#71818B"/>')
def place(k,x,y,slot):
 d=D[k];xx=x+(slot-d.width)/2;yy=y-d.height/2
 parts.append(f'<g transform="translate({xx},{yy})">{d.svg_body}</g>');geo=d.translated_geometry(xx,yy)
 elements.append(FigureElement(kind='molecule',label=k,box=geo['ink_box'],label_boxes=geo['label_boxes'],bond_lengths_px=geo['bond_lengths_px']))
text(width/2,27,'FDAA-assisted chiral derivatization and TIMS-MS discrimination of alanine enantiomers','title')
text(24,58,'FDAA chiral derivatization','heading','start')
ax=(width-(aw+fw+pw+196))/2;fx=ax+aw+48;arrx=fx+fw+18;px=arrx+130
for row,(aa,p,config) in enumerate([('D-Ala','L-FDAA-D-Ala','S,R'),('L-Ala','L-FDAA-L-Ala','S,S')]):
 y=80+rowh/2+row*rowh
 place(aa,ax,y,aw);text(ax+aw+24,y+5,'+')
 place('L-FDAA',fx,y,fw);arrow(arrx,y,arrx+110);text(arrx+55,y-12,'SNAr','small')
 place(p,px,y,pw)
 labely=y+rowh/2-12
 text(ax+aw/2,labely,aa+' ('+('R' if row==0 else 'S')+')','label '+('blue' if row==0 else 'orange'))
 text(fx+fw/2,labely,'L-FDAA (S)','label')
 text(px+pw/2,labely,'L-FDAA–'+aa+' ('+config+')','label '+('blue' if row==0 else 'orange'))
base=80+2*rowh
text(width/2,base+15,'Diastereomers: same neutral formula C₁₂H₁₅N₅O₇ and molecular mass','small')
text(24,base+43,'TIMS-MS discrimination','heading','start')
# Qualitative schematic with deliberately unassigned color-to-peak correspondence.
text(136,base+85,'Diastereomer formation','label');text(136,base+108,'S,R versus S,S','small');arrow(250,base+91,width/2-160)
plotx=width/2-120;plotw=240;plotbase=base+143;ploth=64
parts.append(f'<path class="axis" d="M {plotx} {plotbase-ploth} V {plotbase} H {plotx+plotw}"/>')
for center in (.34,.68):
 pts=[]
 for i in range(241):
  t=i/240;yy=plotbase-ploth*.85*math.exp(-((t-center)/.066)**2/2);pts.append(f'{plotx+t*plotw:.2f},{yy:.2f}')
 parts.append('<polyline points="'+' '.join(pts)+'" stroke="#596C7D" stroke-width="1.8" fill="none"/>')
text(plotx+plotw/2,base+63,'Different ion mobilities','label')
text(plotx+plotw/2,plotbase+20,'Inverse reduced mobility, 1/K₀','small')
parts.append(f'<text x="{plotx-15}" y="{plotbase-ploth/2}" text-anchor="middle" class="small" transform="rotate(-90 {plotx-15} {plotbase-ploth/2})">Relative intensity</text>')
arrow(plotx+plotw+22,base+91,width-255)
text(width-135,base+85,'Enantiomer discrimination','label');text(width-135,base+108,'Via diastereomeric derivatives','small')
text(width/2,base+188,'Schematic illustration; peak order unassigned. TIMS measures derivative ions, potentially metal-adduct ions.','small')
text(width/2,base+207,'Xie et al., Chem. Sci. 2022 • DOI: 10.1039/D2SC03604E; Marfey 1984 • DOI: 10.1007/BF02908688','small')
parts.append('</svg>');root=ET.fromstring('\n'.join(parts));meta=ET.SubElement(root,'{http://www.w3.org/2000/svg}metadata')
meta.text=json.dumps({'structures':S,'validation':{'CIP':expected,'products_formula':'C12H15N5O7','exact_neutral_mass_Da':rdMolDescriptors.CalcExactMolWt(ref),'F_replaced_at_original_aryl_carbon':ar,'original_FDAA_subgraph_preserved':True,'SNAr_balanced_with_formal_HF':True},'provenance':'User supplied stereospecific SMILES and literature citations; paper and PubChem requests HTTP proxy 403; experimental figures not independently accessed','reaction':'Overall SNAr reaction scheme, not elementary curved-arrow mechanism; byproducts omitted graphically, HF included in bookkeeping only','TIMS':'Schematic, no numerical mobility, resolution or peak-order assignments; neutral structures depict reaction products, actual TIMS ions/adduct unspecified','conditions':'Optional literature conditions omitted pending direct verification'},ensure_ascii=False)
ET.register_namespace('','http://www.w3.org/2000/svg');out=Path(__file__).resolve().with_suffix('.svg');out.write_text(ET.tostring(root,encoding='unicode'));png=out.with_suffix('.png')
assert not svg_to_png(out,png,scale=4800/width)
with Image.open(png) as im:im.load();im.save(png,dpi=(600,600))
warn,ink=raster_ink_check(png);assert not warn and ink['status']=='ok'
audit=audit_figure(renderer='rdkit_FDAA_SNAr_TIMS',width=width,height=height,elements=elements,target_bond_px=30,png_path=png);assert not audit['blocking_issues'] and not audit['warnings'],audit
assert out.stat().st_size>0 and png.stat().st_size>0;ET.parse(out)
print('Canvas',width,height,'PNG',ink['size'],'validated',audit['warnings'])

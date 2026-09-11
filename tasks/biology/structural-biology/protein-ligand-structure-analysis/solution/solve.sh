#!/usr/bin/env bash
set -euo pipefail

ROOT="${TASK_ROOT:-/app}"
SEQ="${ROOT}/sequence.fasta"
PDB="${ROOT}/template.pdb"
SDF="${ROOT}/ligand_library.sdf"
OUTDIR="${ROOT}/solution_output"

mkdir -p "$OUTDIR"

python3 - "$SEQ" "$PDB" "$SDF" "$OUTDIR" <<'PY'
import sys, os, re, json, math
from collections import defaultdict

seq_path, pdb_path, sdf_path, outdir = sys.argv[1:]
os.makedirs(outdir, exist_ok=True)

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name=line[1:].split()[0]
                buf=[]
            else:
                buf.append(re.sub(r"[^A-Za-z]", "", line).upper())
    if name is not None: seqs[name]="".join(buf)
    return seqs

def nw(a,b,match=1,mismatch=-1,gap=-2):
    n,m=len(a),len(b)
    dp=[[0]*(m+1) for _ in range(n+1)]
    tb=[[None]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=i*gap; tb[i][0]="U"
    for j in range(1,m+1): dp[0][j]=j*gap; tb[0][j]="L"
    for i in range(1,n+1):
        for j in range(1,m+1):
            vals=[
                (dp[i-1][j-1]+(match if a[i-1]==b[j-1] else mismatch),"D"),
                (dp[i-1][j]+gap,"U"),
                (dp[i][j-1]+gap,"L")]
            dp[i][j],tb[i][j]=max(vals,key=lambda x:x[0])
    aa=[];bb=[];i=n;j=m
    while i or j:
        t=tb[i][j]
        if t=="D": aa.append(a[i-1]); bb.append(b[j-1]); i-=1;j-=1
        elif t=="U": aa.append(a[i-1]); bb.append("-"); i-=1
        else: aa.append("-");bb.append(b[j-1]);j-=1
    return "".join(reversed(aa)),"".join(reversed(bb)),dp[n][m]

def identity(aln_a,aln_b):
    cols=[(x,y) for x,y in zip(aln_a,aln_b) if x!="-" and y!="-"]
    return (sum(x==y for x,y in cols)/len(cols)*100 if cols else 0.0),len(cols)

def parse_pdb(path):
    atoms=[]
    seqres=defaultdict(list)
    het=[]
    missing=[]
    resolution=None; rwork=None; rfree=None; mean_b=None; esu=None
    bond_rms=None; angle_rms=None
    with open(path,errors="replace") as f:
        for line in f:
            rec=line[:6].strip()
            if rec=="SEQRES":
                chain=line[11].strip()
                toks=line[19:70].split()
                seqres[chain].extend(toks)
            elif rec in ("ATOM","HETATM"):
                try:
                    x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
                except: continue
                atom={
                    "record":rec,"atom":line[12:16].strip(),"res":line[17:20].strip(),
                    "chain":line[21].strip(),"resi":line[22:26].strip(),"x":x,"y":y,"z":z,
                    "elem":line[76:78].strip().upper() or re.sub(r"[^A-Z]","",line[12:16]).upper()[:1]
                }
                atoms.append(atom)
                if rec=="HETATM" and atom["res"] not in {"HOH","WAT","DOD"}:
                    het.append(atom)
            elif line.startswith("REMARK   2 RESOLUTION."):
                m=re.search(r"RESOLUTION\.\s+([0-9.]+)",line)
                if m: resolution=float(m.group(1))
            elif "R VALUE" in line and "WORKING" in line and "FREE" in line:
                nums=re.findall(r"0\.\d+",line)
                if len(nums)>=2: rwork=float(nums[0]); rfree=float(nums[1])
            elif "ESTIMATED OVERALL COORDINATE ERROR" in line:
                esu=line.strip()
            elif "RMS BOND DISTANCE" in line:
                m=re.search(r"([0-9.]+)\s+ANGSTROMS",line)
                if m: bond_rms=float(m.group(1))
            elif "RMS BOND ANGLE" in line:
                m=re.search(r"([0-9.]+)\s+DEGREES",line)
                if m: angle_rms=float(m.group(1))
    return atoms,seqres,het,{"resolution":resolution,"rwork":rwork,"rfree":rfree,
                              "coordinate_error":esu,"bond_rms":bond_rms,"angle_rms":angle_rms}

def parse_sdf(path):
    records=[]
    with open(path,errors="replace") as f: text=f.read()
    for block in text.split("$$$$"):
        if not block.strip(): continue
        lines=block.splitlines()
        name=lines[0].strip() if lines else ""
        props={}
        # Parse V2000 atoms and simple properties. Coordinates are retained for frame check.
        nat=0; coords=[]
        if len(lines)>3:
            try:
                nat=int(lines[3][0:3])
            except: nat=0
            for line in lines[4:4+nat]:
                try:
                    coords.append((float(line[0:10]),float(line[10:20]),float(line[20:30]),
                                   line[31:34].strip()))
                except: pass
        i=4+nat
        while i<len(lines):
            if lines[i].startswith(">"):
                m=re.search(r"<([^>]+)>",lines[i])
                key=m.group(1) if m else ""
                i+=1; vals=[]
                while i<len(lines) and lines[i].strip() and not lines[i].startswith(">"):
                    vals.append(lines[i].strip()); i+=1
                props[key]=" ".join(vals)
            else: i+=1
        records.append({"name":name,"props":props,"coords":coords})
    return records

seqs=read_fasta(seq_path)
target=next(iter(seqs.values())) if seqs else ""
atoms,seqres,het,quality=parse_pdb(pdb_path)
template_seq="".join({"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
"HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T",
"TRP":"W","TYR":"Y","VAL":"V"}.get(x,"X") for x in seqres[next(iter(seqres))]) if seqres else ""

aln_t,aln_p,score=nw(target,template_seq)
ident,ncomp=identity(aln_t,aln_p)

# Identify structural copies and ligand-defined pockets.
protein=[a for a in atoms if a["record"]=="ATOM"]
chains=sorted(set(a["chain"] for a in protein))
lig_res=defaultdict(list)
for a in het: lig_res[(a["chain"],a["res"],a["resi"])].append(a)

# For each non-water HET group, collect nearby protein residues at <=4 A.
pockets=[]
for key,latoms in lig_res.items():
    if key[1] in {"HOH","WAT","DOD"}: continue
    contacts={}
    for la in latoms:
        for pa in protein:
            d=math.dist((la["x"],la["y"],la["z"]),(pa["x"],pa["y"],pa["z"]))
            if d<=4.0:
                rk=(pa["chain"],pa["resi"],pa["res"])
                contacts[rk]=min(d,contacts.get(rk,999))
    if contacts:
        pockets.append({"ligand_group":key,"n_contacts":len(contacts),
                        "residues":sorted([{"chain":c,"resi":r,"resname":n,"min_distance_A":round(d,3)}
                                           for (c,r,n),d in contacts.items()],
                                          key=lambda x:x["min_distance_A"])})
pockets.sort(key=lambda x:x["n_contacts"],reverse=True)

# Heme-proximal structural pocket: residues around heme non-hydrogen atoms.
heme=[a for a in het if a["res"]=="HEM"]
heme_contacts={}
if heme:
    for la in heme:
        for pa in protein:
            d=math.dist((la["x"],la["y"],la["z"]),(pa["x"],pa["y"],pa["z"]))
            if d<=6.0:
                rk=(pa["chain"],pa["resi"],pa["res"])
                heme_contacts[rk]=min(d,heme_contacts.get(rk,999))

ligands=parse_sdf(sdf_path)

# Detect whether ligand coordinates plausibly share the protein coordinate frame.
# This is deliberately conservative: absolute centroid proximity alone is not treated as docking.
pcent=(sum(a["x"] for a in protein)/len(protein),sum(a["y"] for a in protein)/len(protein),
       sum(a["z"] for a in protein)/len(protein)) if protein else None
for L in ligands:
    c=None
    if L["coords"]:
        c=(sum(x for x,y,z,e in L["coords"])/len(L["coords"]),
           sum(y for x,y,z,e in L["coords"])/len(L["coords"]),
           sum(z for x,y,z,e in L["coords"])/len(L["coords"]))
    L["coordinate_centroid"]=c
    L["coordinate_frame_comparable"]=bool(c and pcent and math.dist(c,pcent)<80.0)

# Extract key annotation fields without inventing missing values.
def pick(props,*keys):
    for k in keys:
        if k in props and props[k]!="": return props[k]
    return None

lig_summary=[]
for L in ligands:
    p=L["props"]
    lig_summary.append({
        "record_name":L["name"],
        "bindingdb_monomerid":pick(p,"BindingDB MonomerID"),
        "ligand_name":pick(p,"BindingDB Ligand Name"),
        "chembl_id":pick(p,"ChEMBL ID of Ligand"),
        "pubchem_cid":pick(p,"PubChem CID of Ligand"),
        "inchi_key":pick(p,"Ligand InChI Key"),
        "canonical_smiles":pick(p,"Canonical_SMILES"),
        "target_name":pick(p,"Target Name"),
        "target_organism":pick(p,"Target Source Organism According to Curator or DataSource"),
        "ki_nM":pick(p,"Ki (nM)"),
        "ic50_nM":pick(p,"IC50 (nM)"),
        "kd_nM":pick(p,"Kd (nM)"),
        "source":pick(p,"Curation/DataSource"),
        "doi":pick(p,"Article DOI"),
        "n_atoms":len(L["coords"]),
        "coordinate_frame_comparable":L["coordinate_frame_comparable"]
    })

# Ranking is intentionally conservative. If supplied SDF coordinates are not demonstrably
# in the PDB frame, no geometry-derived ligand ranking is produced.
frame_ok=all(x["coordinate_frame_comparable"] for x in ligands) if ligands else False
decision=("NO_DEFENSIBLE_STRUCTURAL_RANKING" if not frame_ok
          else "STRUCTURAL_COMPARISON_POSSIBLE_BUT_REQUIRES_EXPLICIT_PLACEMENT_VALIDATION")

result={
 "sequence_comparison":{
   "target_length":len(target),"template_length":len(template_seq),
   "alignment_method":"global Needleman-Wunsch; match=+1, mismatch=-1, gap=-2",
   "alignment_score":score,"aligned_non_gap_pairs":ncomp,"sequence_identity_percent":round(ident,3),
   "alignment_target":aln_t,"alignment_template":aln_p
 },
 "structure":{
   "chains":chains,"quality":quality,
   "ligand_defined_pockets":pockets[:10],
   "heme_proximal_residue_count":len(heme_contacts),
   "heme_proximal_residues":sorted([{"chain":c,"resi":r,"resname":n,"min_distance_A":round(d,3)}
                                     for (c,r,n),d in heme_contacts.items()],
                                    key=lambda x:x["min_distance_A"])[:100]
 },
 "ligands":lig_summary,
 "contact_criterion":"<=4.0 A between ligand non-hydrogen atoms and protein non-hydrogen atoms",
 "contact_analysis_status":"not_performed_as_binding_contacts" if not frame_ok else "geometry_available_for_validation",
 "interpretation_limitations":[
   "SDF ligand coordinates are not assumed to be superposed into the PDB coordinate frame.",
   "Absence of protein-ligand contacts is therefore not interpreted as incompatibility.",
   "Ligand annotations such as Ki are reported as supplied records and are not used as proof of structural binding.",
   "A template-derived pocket is a structural hypothesis, not experimental proof of affinity or activity."
 ],
 "final_assessment":decision
}

with open(os.path.join(outdir,"results.json"),"w") as f:
    json.dump(result,f,indent=2)

with open(os.path.join(outdir,"report.md"),"w") as f:
    f.write("# Protein–ligand structural investigation\n\n")
    f.write("## Sequence/template comparison\n")
    f.write(f"- Target length: {len(target)} residues\n- Template sequence length represented by SEQRES: {len(template_seq)} residues\n")
    f.write(f"- Global alignment: Needleman–Wunsch (match +1, mismatch −1, gap −2)\n")
    f.write(f"- Sequence identity over aligned non-gap pairs: **{ident:.3f}%** ({ncomp} pairs)\n\n")
    f.write("The template is evaluated as a structural analogue rather than as an identical target model. "
            "Sequence differences and the template's crystallographic uncertainty should be considered when transferring "
            "pocket geometry or residue-level interaction hypotheses.\n\n")
    f.write("## Structural quality and pocket\n")
    for k,v in quality.items():
        f.write(f"- {k}: {v}\n")
    f.write(f"- Protein chains observed in ATOM records: {', '.join(chains) or 'none'}\n")
    if pockets:
        f.write(f"- Heteroatom-defined contact groups with protein contacts at ≤4.0 Å: {len(pockets)}\n")
        for p in pockets[:5]:
            f.write(f"  - {p['ligand_group']}: {p['n_contacts']} protein residues within 4.0 Å\n")
    else:
        f.write("- No non-water heteroatom group yielded a ≤4.0 Å protein contact set.\n")
    f.write("\n")
    f.write("## Ligand library\n")
    for x in lig_summary:
        f.write(f"- **{x['ligand_name'] or x['record_name']}**; ChEMBL={x['chembl_id']}; "
                f"BindingDB={x['bindingdb_monomerid']}; Ki={x['ki_nM']} nM; target={x['target_name']}; "
                f"source={x['source']}\n")
    f.write("\n## Contact analysis limitation\n")
    if not frame_ok:
        f.write("The supplied ligand coordinates are not treated as experimentally placed in the template's coordinate "
                "frame. Consequently, a ≤4.0 Å protein–ligand contact calculation is not used to claim incompatibility "
                "or compatibility. A placement/superposition step would be required before geometric contact ranking.\n\n")
    else:
        f.write("The coordinate sets pass the conservative frame check, but this alone does not establish a binding pose; "
                "contacts should be interpreted as geometric observations only.\n\n")
    f.write("## Final assessment\n")
    f.write(f"**{decision}**\n\n")
    f.write("No compound is promoted solely because it has a favorable database annotation. Structural observations are "
            "kept separate from supplied activity measurements, and neither is treated as proof of biological activity.\n")
PY

echo "Results written to ${OUTDIR}/results.json"
echo "Report written to ${OUTDIR}/report.md"


"""Construct the molecular event roster from the supplementary S1/S2 workbook."""
import re
from pathlib import Path
import pandas as pd


def sheet_records(path):
    from openpyxl import load_workbook
    workbook=load_workbook(path,read_only=True,data_only=True)
    def section(prefix,label):
        values=list(next(s for s in workbook if s.title.startswith(prefix)).values)
        header_index=next(i for i,row in enumerate(values) if row and str(row[0]).startswith(label))+2
        rows=[]
        for row in values[header_index+1:]:
            if not row or not row[0] or (prefix=='S2' and not re.fullmatch(r'SGR[0-9]+',str(row[0]))):break
            rows.append(row)
        return values[header_index],rows
    return [section('S1','A.'),section('S2','B.')]


def roster(path,archive):
    sections=sheet_records(path)
    patients=pd.DataFrame(sections[0][1],columns=sections[0][0]).set_index('Patient ID')
    source=pd.DataFrame(sections[1][1],columns=sections[1][0])
    modeldir=archive/'analysis'
    previous=pd.read_csv(modeldir/'source_data/S24_Structural_pairs.csv').set_index('SGR')
    old_events=pd.read_csv(modeldir/'source_data/S3_Second_hits.csv').set_index('Somatic_Event_ID')
    revised = pd.read_csv(Path(__file__).resolve().parents[2] / 'data/protein_normalizations.csv')
    revised = revised.loc[revised.variant_id.str.startswith('SOMVAR')].set_index('normalized_cdna')
    aa3=dict(zip('ACDEFGHIKLMNPQRSTVWY','Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln Arg Ser Thr Val Trp Tyr'.split()))
    def formatted_germline(value):
        parsed=re.fullmatch(r'([A-Z])(\d+)([A-Z])',value)
        return ('p.'+aa3[parsed[1]]+parsed[2]+aa3[parsed[3]],int(parsed[2])) if parsed else (None,None)
    records=[]
    for _,row in source.iterrows():
        ev=row['Event ID']; patient=row['Patient ID']; bg=row['Germline background']
        phase=row['Phase status']
        phase={'trans_excluded':'confirmed_trans','confirmed_cis_and_trans':'confirmed_cis_both'}.get(phase,phase)
        if ev in previous.index:
            original=previous.loc[ev]
            g_label,g_pos=original['Germline'],original['Germline aa']
            s_label,s_pos=original['Second hit'],original['Second-hit aa']
        else:
            g_label,g_pos=formatted_germline(bg)
            s_label=row['Normalized protein']
            s_pos=int(old_events.loc[ev,'Analysis_position_aa']) if ev in old_events.index else int(re.search(r'\d+',s_label).group())
        if row['Acquired cDNA'] in revised.index:
            normalized = revised.loc[row['Acquired cDNA']]
            s_label, s_pos = normalized.normalized_protein, normalized.residue
        seq=row['Included: fixed-window RNA predictions']=='Yes'
        assert seq==(row['Included: gene-centered RNA / chromatin predictions']=='Yes')
        struct=row['Included: residue-based structural analysis']=='Yes'
        reason=''
        if not seq:
            if phase=='confirmed_trans':reason='demonstrated trans'
            elif phase=='dual_germline_excluded':reason='two germline variants in the same individual'
            elif pd.isna(row['Acquired cDNA']):reason='exact cDNA unknown and combined-sequence prediction unavailable'
            else: reason='combined-sequence prediction unavailable'
        records.append(dict(event=ev,patient=patient,source_patient=patients.loc[patient,'Source patient ID'],germline=bg,somatic_cdna=row['Acquired cDNA'],
                            reported_protein=row['Protein as reported'],germline_label=g_label,somatic_label=s_label,germline_pos=g_pos,somatic_pos=s_pos,phase=phase,
                            method=row['Phase evidence method'],structural=struct,sequence=seq,reason_without_sequence=reason))
    return pd.DataFrame(records)

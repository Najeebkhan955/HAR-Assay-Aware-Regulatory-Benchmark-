from pathlib import Path
import csv
from artifact_tool import Workbook, SpreadsheetFile
OUT=Path('/mnt/data/HAR_Stage16_Comparative_TF');wb=Workbook.create()
files=[('HAR_Readiness','S1_All180_HAR_Sequence_Mapping_Readiness.csv'),('Bait_Fragment_Map','S2_All_Pal_Bait_Fragments_for_Screen.csv'),('Pilot_Catalog','S3_TenHAR_Provisional_Catalog_Matches.csv'),('TF_Expression','S4_NonTargetingControl_NeuralTF_Expression.csv'),('Pilot_PWM','S5_TenHAR_SixTF_PWM_Exploratory_Pilot.csv'),('Source_Register','S6_Independent_Source_Register.csv')]
for nm,fn in files:
 with open(OUT/fn,encoding='utf-8',newline='') as f:
  data=list(csv.reader(f))
 if nm=='Bait_Fragment_Map':
  # the full bait table remains a standalone CSV; workbook contains limited preview to remain responsive
  data=data[:351]
 sht=wb.worksheets.add(nm)
 # normalize numeric and true/false strings for filtering
 vals=[]
 for row in data:
  if row is data[0]:vals.append(row);continue
  tmp=[]
  for idx,c in enumerate(row):
   hname=data[0][idx].lower()
   if c and any(s in hname for s in ("score","fraction","_cpm","_counts","length_bp","overlap_bp","_cells","number_","_start","_end","_hits","delta_","_batch","_intervals")) and ("status" not in hname):
    try:
     tmp.append(float(c) if "." in c or "e" in c.lower() else int(c));continue
    except ValueError:pass
   if c in ('True','False'):tmp.append(c=='True')
   elif c and c.strip()==c and c.lstrip('-').isdigit() and not c.startswith('0'):
    tmp.append(int(c))
   else:
    try:
     if c and ('.' in c) and not c.startswith('0') and 'http' not in c:tmp.append(float(c))
     else:tmp.append(c)
    except:tmp.append(c)
  vals.append(tmp)
 if not vals:continue
 sht.get_range_by_indexes(0,0,len(vals),len(vals[0])).values=vals
 sht.freeze_panes.freeze_rows(1)
 headers=sht.get_range_by_indexes(0,0,1,len(vals[0]))
 headers.format={'fill':'#17465A','font':{'bold':True,'color':'#FFFFFF'},'row_height':38,'wrap_text':True,'vertical_alignment':'center'}
 for idx,h in enumerate(data[0]):
  col=sht.get_range_by_indexes(0,idx,len(vals),1)
  col.format.column_width=32 if any(k in h.lower() for k in ('interpret','source','limit','evidence','status','meaning','url','catalog','alias')) else (23 if len(h)>18 else 17)
 try:
  sht.tables.add(sht.get_range_by_indexes(0,0,len(vals),len(vals[0])),True,'Table_'+nm.replace('_',''))
 except Exception as e:print('table exception',nm,str(e)[:120])
 print(nm,len(vals),len(vals[0]))
SpreadsheetFile.export_xlsx(wb).save(str(OUT/'HAR_Stage16_Comparative_Sequence_PWM_Supplementary.xlsx'))
print(wb.inspect({'kind':'table','range':'Pilot_PWM!A1:H4','include':'values','table_max_rows':5}).ndjson[:850])

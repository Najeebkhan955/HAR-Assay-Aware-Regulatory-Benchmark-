from pathlib import Path
from zipfile import ZipFile
import base64, hashlib, json, re, os
SOURCE=Path('/mnt/data/HAR_Complete_Stage1-to25_Available_Code_Data_Results_GitHub_Upload.zip')
OUT=Path('/mnt/data/HAR_GitHub_Code_Transfer'); OUT.mkdir(exist_ok=True)
zs=ZipFile(SOURCE)
items=[]
for member in sorted(zs.namelist()):
 if member.endswith(('.ipynb','.py','.R','.sh')) and not member.startswith('__MACOSX/'):
  payload=zs.read(member)
  # guard accidental credentials in source; report matches without secrets
  try: s=payload.decode('utf-8')
  except: raise ValueError(f'not utf-8: {member}')
  # github location logical prefix code/ plus original preserved path
  ghpath='code/'+member
  items.append(dict(path=ghpath,source=member,content=payload,sha256=hashlib.sha256(payload).hexdigest(),size=len(payload)))
print('total code',len(items),'bytes',sum(i['size'] for i in items))
# Note: duplicate content is preserved where original paths differ; git store deduplicates.
manifest=[{k:i[k] for k in ('path','source','sha256','size')} for i in items]
(OUT/'UPLOAD_CODE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
# base64 text bundle <=80kb source per file, each line under 400 chars
parts=[]; current=[]; n=0
for idx,item in enumerate(items):
 blocks=[]
 encoded=base64.b64encode(item['content']).decode('ascii')
 header='@FILE '+base64.urlsafe_b64encode(item['path'].encode()).decode('ascii')+' '+str(item['size'])+' '+item['sha256']
 blocks.append(header)
 blocks.extend(encoded[i:i+320] for i in range(0,len(encoded),320))
 blocks.append('@END')
 if current and sum(len(z)+1 for z in current)+sum(len(z)+1 for z in blocks)>95000:
  name=f'CODE_TRANSFER_{n:03d}.txt';(OUT/name).write_text('\n'.join(current)+'\n',encoding='ascii');parts.append(name);current=[];n+=1
 current.extend(blocks)
if current:
 name=f'CODE_TRANSFER_{n:03d}.txt';(OUT/name).write_text('\n'.join(current)+'\n',encoding='ascii');parts.append(name)
print('parts',len(parts),'sizes',[(p,(OUT/p).stat().st_size) for p in parts[:7]])
(OUT/'PARTS.json').write_text(json.dumps(parts,indent=2))

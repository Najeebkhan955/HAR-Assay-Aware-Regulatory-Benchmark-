"""Write a relative SHA256SUMS.txt, excluding the manifest itself and local cache files."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
exclude={'SHA256SUMS.txt'}
paths=sorted(p for p in root.rglob('*') if p.is_file() and p.name not in exclude and '__pycache__' not in p.parts and 'HAR_Stage23_WTc_N2_Cache' not in p.parts and 'HAR_Stage24_WTc_GEO_Cache' not in p.parts and 'HAR_Stage24_WTc_Results' not in p.parts and 'HAR_Stage24_Embedded' not in p.parts)
rows=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(root).as_posix() for p in paths]
(root/'SHA256SUMS.txt').write_text('\n'.join(rows)+'\n')
print('Wrote',len(rows),'file checksums')

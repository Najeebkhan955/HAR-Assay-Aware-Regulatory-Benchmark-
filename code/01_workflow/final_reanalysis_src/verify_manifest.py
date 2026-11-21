"""Verify all packaged research bytes with SHA256SUMS.txt."""
from pathlib import Path
import hashlib,sys
root=Path(__file__).resolve().parents[1];lines=(root/'SHA256SUMS.txt').read_text().splitlines()
errors=[]
for s in lines:
    dig,rel=s.split('  ',1);p=root/rel
    if not p.is_file():errors.append('MISSING '+rel)
    elif hashlib.sha256(p.read_bytes()).hexdigest()!=dig:errors.append('HASH MISMATCH '+rel)
if errors:
    print('ERRORS:',*errors,sep='\n');sys.exit(1)
print('SHA256 checksum PASS:',len(lines),'files')

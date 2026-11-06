"""Check the integrity of the distributed inputs and code."""
from pathlib import Path
import csv,hashlib
r=Path(__file__).resolve().parent
with (r/'MANIFEST_SHA256.csv').open(newline='') as f:
 rows=list(csv.DictReader(f))
for row in rows:
 p=r/row['path']
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
  raise SystemExit('Integrity check failed: '+row['path'])
print('PASS:',len(rows),'files')

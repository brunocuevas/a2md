import json
from pathlib import Path

path2files = Path(__file__).parent.parent

with open(Path(__file__).parent.parent / 'a2mdtests/gdb_test_files.json') as f:
	gdb_test = json.load(f)

for k, item in gdb_test.items():
	item['path'] = path2files / item['path']

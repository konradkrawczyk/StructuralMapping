import sys

for line in sys.stdin:
  if line.startswith('#'):
    print line.strip()
  else:
    print line[31:38].strip(), line[38:46].strip(), line[46:54].strip()

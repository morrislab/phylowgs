import sys

def redoIDs(fin):
	f = open(fin)
	d = f.readlines()
	f.close()
	head = d[0]
	d = d[1:]
	d = [x.split('\t') for x in d]
	d = [x[1:] for x in d]
	d = [['s%d' % d.index(x)] + x for x in d]
	d = ['\t'.join(x) for x in d]
	f = open(fin,'w')
	f.write(''.join([head]+d))
	f.close()

if __name__ == '__main__':
	redoIDs(sys.argv[1])
	

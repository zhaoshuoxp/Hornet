# Run an original Hornet script with random.seed and a py3 shim for AlignmentFile.next()
import sys, random, runpy, pysam
seed, script = int(sys.argv[1]), sys.argv[2]
class P:
    def __init__(self, *a, **k):
        if isinstance(k.get('template'), P): k['template'] = k['template']._f
        self._f = pysam.AlignmentFile(*a, **k)
    def next(self): return next(self._f)
    def __iter__(self): return iter(self._f)
    def __getattr__(self, n): return getattr(self._f, n)
pysam.Samfile = P
random.seed(seed)
sys.argv = [script] + sys.argv[3:]
runpy.run_path(script, run_name='__main__')

import sys
try: 
    from calign import aligner
    print >>sys.stderr, "USING CYTHON"
except ImportError:
    print >>sys.stderr, "USING PYTHON"
    import matrix
    from align import aligner

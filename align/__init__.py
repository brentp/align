import sys
try: 
    from calign import aligner
except ImportError:
    import matrix
    from align import aligner

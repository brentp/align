try: 
    from calign import aligner, score_alignment
except ImportError:
    from . import matrix
    from .align import aligner

try: 
    from calign import aligner, score_alignment
except ImportError:
    import matrix
    from align import aligner

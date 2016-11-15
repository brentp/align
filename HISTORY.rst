.. :changelog:

History
=======

Version 0.0.3
-------------

Backwards incompatible changes since version 0.0.2:

    * Alignment results are now stored in an object called ``AlignmentResult``.
      The pure Python implementation is based on the ``namedtuple`` class.

    * The ``AlignmentResult`` object, in addition to containing the aligned
      sequences, also contains the alignment start positions, end positions,
      number of gaps on each sequence, number of mismatches, and the alignment
      score.

    * A new argument, ``max_hits``, is added to the ``aligner``. This
      determines the maximum number of alignments to return in case there
      are multiple optimal alignments (i.e. alignment with the same maximum
      score). The default value is 1. When set to ``None``, all optimal
      alignments are returned, except for when the method is ``global``.
      In this case, only one optimal alignment is returned. This also changes
      the return type of ``aligner`` to be a list, regardless of how many
      alignments are returned.

Additionally, several internal code refactor were done.

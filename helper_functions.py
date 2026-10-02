from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

"""Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
def global_alignment(seq1, seq2, scoring_function):
    # Gap penalty
    gap = -1

    n = len(seq1)
    m = len(seq2)

    # ---------------------------------------------------------
    # 1. Create and initialise the dynamic programming matrix
    # ---------------------------------------------------------

    # score[i][j] = best alignment score for
    # seq1[:i] and seq2[:j]
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Aligning a sequence with an empty sequence requires
    # inserting a gap for every character.
    for i in range(1, n + 1):
        score[i][0] = i * gap

    for j in range(1, m + 1):
        score[0][j] = j * gap

    # ---------------------------------------------------------
    # 2. Fill the dynamic programming matrix
    # ---------------------------------------------------------

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            # Option 1: align the two characters
            diagonal = (
                score[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            # Option 2: align seq1[i-1] with a gap
            up = score[i - 1][j] + gap

            # Option 3: align seq2[j-1] with a gap
            left = score[i][j - 1] + gap

            # Keep the best possible score
            score[i][j] = max(diagonal, up, left)

    # ---------------------------------------------------------
    # 3. Traceback
    # ---------------------------------------------------------

    aligned_seq1 = []
    aligned_seq2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        # Diagonal: seq1[i-1] aligned with seq2[j-1]
        if (
            i > 0
            and j > 0
            and score[i][j]
            == score[i - 1][j - 1]
            + scoring_function(seq1[i - 1], seq2[j - 1])
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1

        # Up: seq1[i-1] aligned with a gap
        elif (
            i > 0
            and score[i][j] == score[i - 1][j] + gap
        ):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1

        # Left: seq2[j-1] aligned with a gap
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])

            j -= 1

    # Traceback constructs the sequences backwards,
    # so reverse them.
    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(score[n][m])

"""Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
def local_alignment(seq1, seq2, scoring_function):
    gap = -1
    n = len(seq1)
    m = len(seq2)

    # ---------------------------------------------------------
    # 1. Initialise: first row/column stay 0 for local alignment
    # ---------------------------------------------------------
    score = [[0] * (m + 1) for _ in range(n + 1)]
    best_score = 0
    best_pos = (0, 0)

    # ---------------------------------------------------------
    # 2. Fill the matrix
    # ---------------------------------------------------------
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            # Scores can never drop below 0 (start a new local alignment)
            score[i][j] = max(0, diagonal, up, left)

            if score[i][j] > best_score:
                best_score = score[i][j]
                best_pos = (i, j)

    # ---------------------------------------------------------
    # 3. Traceback from the best cell until the score hits 0
    # ---------------------------------------------------------
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = best_pos

    while score[i][j] > 0:
        if score[i][j] == score[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1]):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(best_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

# generate scoring function based on blosum62 scoring matrix
# set default gap to -1
def blosum62_matrix(x, y, default=-1):
    try:
        result = blosum62[x][y]
        return result
    except KeyError:
        return default

def percent_identity(seq1, seq2):
    score = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 == base2:
            score = score + 1
    return score / len(seq1)
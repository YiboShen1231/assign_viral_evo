def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    import math

    valid_sites = 0
    differences = 0

    for a, b in zip(reference_sequence.upper(), distant_sequence.upper()):
        if a not in "ACGT" or b not in "ACGT":
            continue

        valid_sites += 1
        if a != b:
            differences += 1

    if valid_sites == 0:
        return 0.0

    p = differences / valid_sites
    distance = -0.75 * math.log(1 - 4 * p / 3)

    return distance * valid_sites



def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """

    import math
    valid_sites = 0
    transitions = 0
    transversions = 0

    for a, b in zip(reference_sequence.upper(), distant_sequence.upper()):
        if a not in "ACGT" or b not in "ACGT":
            continue
        valid_sites += 1
        if a != b:
            if (a in "AG" and b in "AG") or (a in "CT" and b in "CT"):
                transitions += 1
            else:
                transversions += 1

    if valid_sites == 0:
        return 0.0

    p = transitions / valid_sites
    q = transversions / valid_sites
    distance = -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q)

    return distance * valid_sites

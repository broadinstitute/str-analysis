import collections
import numpy as np
import random

from rapidfuzz.distance import Hamming, Levenshtein

from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.fasta_utils import normalize_chrom_using_pyfaidx_fasta
from str_analysis.utils.trf_runner import TRFRunner

HAMMING_DISTANCE_METRIC = "hamming"
EDIT_DISTANCE_METRIC = "edit"
DEFAULT_DISTANCE_METRIC = HAMMING_DISTANCE_METRIC

def compute_repeat_purity(
        nucleotide_sequence,
        motif,
        include_partial_repeats=False,
        distance_metric=DEFAULT_DISTANCE_METRIC):
    """This method generates a sequence of pure repeats of the given motif (with length equal to the length of the input
    nucleotide sequence), and then computes the number of substitutions (ie. interruptions) in the given nucleotide
    sequence, compared to the synthetic pure repeat sequence. If the nucleotide sequence length is not an
    exact multiple of the motif length, the sequence will be trimmed to the nearest multiple of the motif length before
    computing stats, unless include_partial_repeats is True.

    Args:
        nucleotide_sequence (str): nucleotide sequence with A, C, G, T
        motif (str): a repeat motif
        include_partial_repeats (bool): whether to include any partial repeat at the end of the sequence in the
            calculation. If False, only complete repeats will be included.
        distance_metric (HAMMING_DISTANCE_METRIC or EDIT_DISTANCE_METRIC): distance metric to use

    Return:
        2-tuple:
            fraction_pure_bases (float): When HAMMING_DISTANCE_METRIC is used, this fraction represents the
                fraction of bases in the input sequence that matched a same-length pure repeat sequence of the given
                motif. In other words, the number of matches / N  where N is the length of the input sequence.
                If EDIT_DISTANCE_METRIC is used, the fraction represents (N - edit_distance) / N with
                insertions, deletions, and substitutions all adding 1 to the edit distance.
            edit_count (int) The number of substitutions or edits needed to convert the input sequence into a pure
                repeat sequence of the given motif.
    """
    if len(nucleotide_sequence) < len(motif):
        return float('nan'), None

    num_repeats = int(len(nucleotide_sequence)/len(motif))
    if include_partial_repeats:
        pure_sequence = motif.upper() * (num_repeats + 1)
        pure_sequence = pure_sequence[:len(nucleotide_sequence)]
    else:
        pure_sequence = motif.upper() * num_repeats
        nucleotide_sequence = nucleotide_sequence[:len(pure_sequence)]

    if distance_metric == HAMMING_DISTANCE_METRIC:
        distance = Hamming.distance(pure_sequence, nucleotide_sequence)
    elif distance_metric == EDIT_DISTANCE_METRIC:
        # https://rapidfuzz.github.io/Levenshtein/levenshtein.html#distance
        distance = Levenshtein.distance(pure_sequence, nucleotide_sequence, weights=(1, 1, 1))
    else:
        raise ValueError("Unknown distance metric {}".format(distance_metric))

    if distance > len(nucleotide_sequence):
        raise Exception(f"ERROR: edit distance ({distance:,d}) is larger than the nucleotide sequence "
                        f"({len(nucleotide_sequence):,d}bp). Distance metric: {distance_metric}, "
                        f"include_partial_repeats: {include_partial_repeats}, Motif: {motif}\n  {nucleotide_sequence} ")
    fraction_pure_bases = (len(nucleotide_sequence) - distance) / len(nucleotide_sequence)

    return fraction_pure_bases, distance


def compute_motif_purity_for_interval(
        reference_fasta,
        chrom,
        start_0based,
        end_1based,
        motif,
        distance_metric=DEFAULT_DISTANCE_METRIC):
    
    if reference_fasta is None:
        return float('nan'), None

    chrom = f"chr{chrom.replace('chr', '')}"  # make sure chrom has "chr" prefix
    reference_sequence = reference_fasta[chrom][start_0based:end_1based]

    return compute_repeat_purity(reference_sequence, motif, include_partial_repeats=True, distance_metric=distance_metric)


def compute_motif_length_quality(motif_length, motif_length_vs_motif_and_purity):
    optimal_motif_length_purities = []
    other_motif_length_purities = []
    for current_motif_length, (current_motif, purity) in motif_length_vs_motif_and_purity.items():
        if current_motif_length >= motif_length and current_motif_length % motif_length == 0:
            optimal_motif_length_purities.append(purity)
        else:
            other_motif_length_purities.append(purity)

    optimal_motif_length_mean_purity = sum(optimal_motif_length_purities) / len(optimal_motif_length_purities)
    if len(other_motif_length_purities) > 0:
        other_motif_lengths_mean_purity = sum(other_motif_length_purities) / len(other_motif_length_purities)
    else:
        # this is the base-line quality score upper-bound (the average of the distribution is 0.25 and the upper-bound appears to be ~0.3)
        # see the output of
        other_motif_lengths_mean_purity = 0.3

    quality_score = optimal_motif_length_mean_purity - other_motif_lengths_mean_purity

    return max(0, quality_score)


def find_optimal_motif_length(nucleotide_sequence, max_motif_length, distance_metric=DEFAULT_DISTANCE_METRIC, verbose=False):
    """Scan different motif lengths from 1 to max_motif_length to find the one that produces the highest
    repeat purity with respect to the reference sequence at the given locus. For each motif length,
    this method finds the most frequent motif of that length within the reference sequence, then constructs
    a synthetic perfect repeat sequence of that motif (making it the same length as the reference sequence),
    then computes purity of that motif length as the fraction of bases in the reference sequence that match the
    previously constructed perfect repeat sequence of that motif.
    """

    #if verbose:
    #    print("--------------------------------")
    #    print(f"Sequence: {nucleotide_sequence}")
    if not nucleotide_sequence:
        raise ValueError("ERROR: nucleotide_sequence is empty")

    if max_motif_length < 1:
        raise ValueError(f"ERROR: max_motif_length is {max_motif_length}")

    # Compute the purity of each motif length
    motif_length_vs_motif_and_purity = {}
    max_purity = 0
    for motif_length in range(1, max_motif_length+1):
        if motif_length > len(nucleotide_sequence) // 2 and motif_length != len(nucleotide_sequence):
            continue

        motif_length_purity, most_common_motif = compute_motif_length_purity(nucleotide_sequence, motif_length, distance_metric=distance_metric)
        if verbose:
            print(f"{motif_length:3d}bp   purity:  {motif_length_purity:.2f}   {most_common_motif}")

        if most_common_motif is None:
            continue

        motif_length_vs_motif_and_purity[len(most_common_motif)] = (most_common_motif, motif_length_purity)
        max_purity = max(max_purity, motif_length_purity)

    if not motif_length_vs_motif_and_purity:
        return None, float('nan'), float('nan')

    if distance_metric == HAMMING_DISTANCE_METRIC:
        # Compute the quality score of each motif length using the motif_length_vs_motif_and_purity dictionary
        motif_length_vs_motif_and_quality_score = {}
        for motif_length, (motif, purity) in motif_length_vs_motif_and_purity.items():
            if motif_length == 1 and purity < max_purity:
                quality = 0   # avoid an edge case where motif quality is unreasonably high for low-purity homopolymers
            else:
                quality = compute_motif_length_quality(motif_length, motif_length_vs_motif_and_purity)
            motif_length_vs_motif_and_quality_score[motif_length] = (motif, quality)
            if verbose:
                print(f"{motif_length:3d}bp   purity:  {purity:.2f}    quality: {quality}")

        # Find the motif length with the highest quality score
        optimal_motif_length = max(
            motif_length_vs_motif_and_quality_score, key=lambda motif_length: (motif_length_vs_motif_and_quality_score[motif_length][1], -motif_length))
    elif distance_metric == EDIT_DISTANCE_METRIC:
        optimal_motif_length = max(
            motif_length_vs_motif_and_purity, key=lambda motif_length: (motif_length_vs_motif_and_purity[motif_length][1], -motif_length))

        # quality scores don't work well with edit distance, so just select based on purity (aka. minimal edit distance)
    else:
        raise ValueError(f"Unknown distance metric {distance_metric}")

    optimal_motif, optimal_purity = motif_length_vs_motif_and_purity[optimal_motif_length]
    optimal_motif_length_quality_score = compute_motif_length_quality(optimal_motif_length, motif_length_vs_motif_and_purity)

    simplified_optimal_motif, _, _ = find_repeat_unit_without_allowing_interruptions(optimal_motif, allow_partial_repeats=False)
    #if simplified_optimal_motif != optimal_motif:
    #    print(f"WARNING: Simplified optimal motif {simplified_optimal_motif} != optimal motif {optimal_motif} for sequence: {nucleotide_sequence}")  # this happens occasionally due to edge cases

    if verbose:
        null_quality = compute_motif_null_quality_score_for_sequence_length(len(nucleotide_sequence), distance_metric=distance_metric)
        print(f"Optimal motif: {len(optimal_motif)}bp   purity: {optimal_purity:.2f}   quality: {optimal_motif_length_quality_score}   (null quality is: {null_quality}). Scanned "
              f"{len(motif_length_vs_motif_and_purity):,d} motif lengths, their average purity was: "
              f"{sum([x[1] for x in motif_length_vs_motif_and_purity.values()]) / len(motif_length_vs_motif_and_purity):.2f}")

    return simplified_optimal_motif, optimal_purity, optimal_motif_length_quality_score


def find_optimal_motif_length_for_interval(
        pyfaidx_reference_fasta_obj,
        chrom,
        start_0based,
        end_1based,
        max_motif_length,
        distance_metric=DEFAULT_DISTANCE_METRIC,
        verbose=False):

    if pyfaidx_reference_fasta_obj is None:
        return None, float('nan'), float('nan')

    chrom = normalize_chrom_using_pyfaidx_fasta(pyfaidx_reference_fasta_obj, chrom)
    reference_sequence = pyfaidx_reference_fasta_obj[chrom][start_0based:end_1based]

    return find_optimal_motif_length(
        nucleotide_sequence=reference_sequence,
        max_motif_length=min(max_motif_length, end_1based - start_0based),
        distance_metric=distance_metric,
        verbose=verbose)


def compute_motif_length_purity(nucleotide_sequence, motif_length, distance_metric=DEFAULT_DISTANCE_METRIC):
    """Find the most frequent motif of the given length in the given nucleotide sequence, then return the
    purity of the input nucleotide sequence, along with the motif itself.
    """

    # slice the reference sequence into subsequences of length motif_length and then get the most common motif
    if motif_length > len(nucleotide_sequence):
        return float('nan'), None
    elif motif_length == len(nucleotide_sequence):
        return 1, nucleotide_sequence
    else:
        end_index = len(nucleotide_sequence) - len(nucleotide_sequence) % motif_length

    sliced_motif_list = [nucleotide_sequence[i:i+motif_length] for i in range(0, end_index, motif_length)]
    if len(sliced_motif_list) == 0:
        return float('nan'), None

    most_common_motif = collections.Counter(sliced_motif_list).most_common(1)[0][0]

    # remove the first occurrence of the most common motif from the sliced nucleotide sequence list. If this wasn't
    # done, the purity would be artificially skewed higher for longer motif lengths since a larger part of the
    # input sequence would automatically be an exact match with the selected motif sequence - eg. if the motif was
    # the length of half the input sequence, this function would take the 1st or 2nd half of the input sequence
    # as the motif, and then purity would be at least 50% regardless of what was in the rest of the sequence.
    sliced_motif_list.remove(most_common_motif)
    if len(sliced_motif_list) == 0:
        return float('nan'), None

    remaining_nucleotide_sequence = "".join(sliced_motif_list)
    most_common_motif_purity, _ = compute_repeat_purity(
        remaining_nucleotide_sequence, most_common_motif, include_partial_repeats=True, distance_metric=distance_metric)

    return most_common_motif_purity, most_common_motif


def compute_motif_length_purity_for_interval(pyfaidx_reference_fasta_obj, chrom, start_0based, end_1based, motif_length, distance_metric=DEFAULT_DISTANCE_METRIC):
    if pyfaidx_reference_fasta_obj is None:
        return float('nan'), None

    chrom = normalize_chrom_using_pyfaidx_fasta(pyfaidx_reference_fasta_obj, chrom)
    reference_sequence = pyfaidx_reference_fasta_obj[chrom][start_0based:end_1based]

    return compute_motif_length_purity(reference_sequence, motif_length, distance_metric=distance_metric)


def generate_motif_null_distributions_using_random_sequences(min_sequence_length=9, max_sequence_length=1000, distance_metric=DEFAULT_DISTANCE_METRIC, verbose=False):
    random.seed(1)
    sequence_length_to_max_purity_and_quality_score = collections.defaultdict(int)
    current_sequence_length = min_sequence_length
    while current_sequence_length < max_sequence_length:
        max_purity = 0
        max_quality = 0

        if current_sequence_length < 100:
            n_trials = 1000
        elif current_sequence_length < 500:
            n_trials = 200
        elif current_sequence_length < 2000:
            n_trials = 100
        else:
            n_trials = 50

        for _ in range(n_trials): # N trials for each length
            random_sequence = ''.join(random.choices("ACGT", k=current_sequence_length))
            optimal_motif, motif_purity, quality = find_optimal_motif_length(random_sequence, current_sequence_length//2, distance_metric=distance_metric, verbose=False)
            max_purity = max(max_purity, motif_purity)
            max_quality = max(max_quality, quality)

        if verbose:
            print(f"{current_sequence_length:4d}bp   Max purity: {max_purity:.2f}     Max quality score: {max_quality:.2f}")

        sequence_length_to_max_purity_and_quality_score[current_sequence_length] = (max_purity, max_quality)

        if current_sequence_length < 100:
            current_sequence_length += 1
        elif current_sequence_length < 200:
            current_sequence_length += 20
        elif current_sequence_length < 500:
            current_sequence_length += 50
        elif current_sequence_length < 1000:
            current_sequence_length += 100
        else:
            current_sequence_length += 1000


    return sequence_length_to_max_purity_and_quality_score


def compute_motif_null_quality_score_for_sequence_length(sequence_length, distance_metric=DEFAULT_DISTANCE_METRIC):
    """If you generate a random nucleotide sequence of the given length and then run 
    find_optimal_motif_length(random_sequence, max_motif_length = len(random_sequence)/2), 
    the maximum possible motif quality you would expect to get is the number returned by this function.
    This can be used to set a minimum threshold for the quality of the motif length returned by 
    find_optimal_motif_length(..) for a given sequence. If the quality is below (or close to) that, it's essentially random.  
    
    The function uses a best-fit power law function to produce the value very close to the empirical distribution 
    generated by running generate_motif_null_distributions_using_random_sequences(max_sequence_length=10001). 
    The fitted curve matches the empirical function table with an R^2 correlation = 0.96
    """

    if distance_metric == HAMMING_DISTANCE_METRIC:
        a = 3.34850673
        b = 0.52020827
        c = -0.00196589
        k = 7.05754591
        null_quality_score = a * (np.asarray(sequence_length) + k)**(-b) + c    # use a scalar numpy array
    elif distance_metric == EDIT_DISTANCE_METRIC:
        a = 22413.5077
        b = 2.79524050
        c = 0.178207196
        k = 38.5197252
        null_quality_score = a * (np.asarray(sequence_length) + k)**(-b) + c    # use a scalar numpy array
    else:
        raise ValueError(f"Unknown distance metric: {distance_metric}")

    return null_quality_score


def find_optimal_motif_using_TRF(trf_executable_path, nucleotide_sequence, max_motif_length, verbose=False):
    """Runs TandemRepeatFinder on the given sequence and returns a motif if it spans the entire sequence
    (allowing 1 partial repeat at either end), or None if no such repeat is found.
    """
    if nucleotide_sequence < 12:
        return None, float('nan'), float('nan')

    max_score = len(nucleotide_sequence) * 2
    optimal_motif = None
    optimal_motif_score = None

    trf_runner = TRFRunner(trf_executable_path,
                           html_mode=True,
                           match_score = 2,
                           mismatch_penalty = 7,
                           indel_penalty = 7,
                           pm = 80,
                           pi = 10,
                           minscore = 24,
                           min_motif_size=1,
                           max_motif_size=max_motif_length)

    for result in trf_runner.run_TRF_on_nucleotide_sequence(nucleotide_sequence):
        motif_length = len(result["repeat_unit"])
        if result["start_0based"] < motif_length and result["end_1based"] > len(nucleotide_sequence) - motif_length and result["repeat_count"] > 1:
            if optimal_motif_score is None or result["alignment_score"] > optimal_motif_score:
                optimal_motif = result["repeat_unit"]
                optimal_motif_score = result["alignment_score"]
                if optimal_motif_score == max_score:
                    break
                if optimal_motif_score > max_score:
                    raise Exception(f"TRF score {motif_length} is greater than max score {max_score} for sequence {nucleotide_sequence}")

    if optimal_motif is None:
        return None, float('nan'), float('nan')

    fraction_pure_bases, _ = compute_repeat_purity(nucleotide_sequence, optimal_motif, include_partial_repeats=True)
    quality_score = optimal_motif_score / max_score

    return optimal_motif, fraction_pure_bases, quality_score

"""
generate_motif_null_distributions_using_random_sequences(max_sequence_length=10001, distance_metric=HAMMING_DISTANCE_METRIC, verbose=True)
generate_motif_null_distributions_using_random_sequences(max_sequence_length=10001, distance_metric=EDIT_DISTANCE_METRIC, verbose=True)

Null distributions (ie. observed upper-bound values in simulated random sequences of the given length)
This table is stored in str-analysis/str_analysis/data/optimal_motif_null_distribution.json

sequence_length hamming_distance_null_purity edit_distance_null_purity hamming_distance_null_quality edit_distance_null_quality
            9bp                         1.00                      1.00                          0.65                       0.60
           10bp                         1.00                      1.00                          0.84                       0.61
           11bp                         1.00                      1.00                          0.86                       0.68
           12bp                         1.00                      1.00                          0.74                       0.57
           13bp                         0.89                      0.83                          0.65                       0.54
           14bp                         1.00                      1.00                          0.77                       0.55
           15bp                         0.86                      0.86                          0.68                       0.44
           16bp                         0.86                      0.86                          0.70                       0.51
           17bp                         0.83                      0.83                          0.58                       0.43
           18bp                         0.86                      0.86                          0.58                       0.47
           19bp                         0.86                      0.86                          0.61                       0.46
           20bp                         0.86                      0.86                          0.58                       0.47
           21bp                         0.75                      0.75                          0.53                       0.37
           22bp                         0.78                      0.75                          0.56                       0.36
           23bp                         0.88                      0.88                          0.63                       0.47
           24bp                         0.80                      0.80                          0.48                       0.37
           25bp                         0.80                      0.80                          0.54                       0.37
           26bp                         0.80                      0.80                          0.50                       0.38
           27bp                         0.80                      0.80                          0.58                       0.36
           28bp                         0.75                      0.73                          0.50                       0.33
           29bp                         0.77                      0.77                          0.52                       0.37
           30bp                         0.73                      0.73                          0.46                       0.31
           31bp                         0.77                      0.79                          0.54                       0.38
           32bp                         0.73                      0.75                          0.53                       0.33
           33bp                         0.67                      0.71                          0.41                       0.33
           34bp                         0.80                      0.80                          0.53                       0.34
           35bp                         0.69                      0.71                          0.47                       0.34
           36bp                         0.77                      0.77                          0.53                       0.31
           37bp                         0.71                      0.71                          0.47                       0.33
           38bp                         0.73                      0.73                          0.50                       0.31
           39bp                         0.67                      0.71                          0.45                       0.31
           40bp                         0.65                      0.70                          0.39                       0.30
           41bp                         0.71                      0.71                          0.46                       0.32
           42bp                         0.70                      0.70                          0.43                       0.24
           43bp                         0.67                      0.71                          0.45                       0.25
           44bp                         0.69                      0.69                          0.44                       0.25
           45bp                         0.65                      0.67                          0.38                       0.24
           46bp                         0.65                      0.68                          0.41                       0.24
           47bp                         0.60                      0.67                          0.35                       0.26
           48bp                         0.65                      0.68                          0.42                       0.29
           49bp                         0.68                      0.72                          0.47                       0.27
           50bp                         0.65                      0.70                          0.41                       0.27
           51bp                         0.67                      0.67                          0.41                       0.26
           52bp                         0.64                      0.65                          0.37                       0.26
           53bp                         0.65                      0.68                          0.38                       0.25
           54bp                         0.64                      0.65                          0.38                       0.25
           55bp                         0.61                      0.74                          0.35                       0.26
           56bp                         0.62                      0.72                          0.37                       0.25
           57bp                         0.67                      0.71                          0.42                       0.25
           58bp                         0.65                      0.67                          0.39                       0.24
           59bp                         0.60                      0.67                          0.36                       0.27
           60bp                         0.63                      0.67                          0.38                       0.29
           61bp                         0.66                      0.67                          0.40                       0.25
           62bp                         0.65                      0.67                          0.38                       0.23
           63bp                         0.57                      0.65                          0.35                       0.22
           64bp                         0.64                      0.70                          0.36                       0.24
           65bp                         0.64                      0.64                          0.40                       0.25
           66bp                         0.58                      0.66                          0.34                       0.21
           67bp                         0.57                      0.65                          0.34                       0.20
           68bp                         0.61                      0.67                          0.36                       0.25
           69bp                         0.58                      0.66                          0.33                       0.23
           70bp                         0.67                      0.67                          0.41                       0.23
           71bp                         0.56                      0.65                          0.31                       0.20
           72bp                         0.59                      0.63                          0.36                       0.21
           73bp                         0.59                      0.64                          0.35                       0.26
           74bp                         0.58                      0.64                          0.32                       0.20
           75bp                         0.62                      0.62                          0.38                       0.19
           76bp                         0.59                      0.67                          0.36                       0.21
           77bp                         0.61                      0.62                          0.33                       0.19
           78bp                         0.61                      0.59                          0.35                       0.21
           79bp                         0.56                      0.64                          0.31                       0.23
           80bp                         0.61                      0.63                          0.35                       0.22
           81bp                         0.57                      0.68                          0.34                       0.22
           82bp                         0.55                      0.66                          0.29                       0.24
           83bp                         0.58                      0.65                          0.32                       0.20
           84bp                         0.59                      0.68                          0.34                       0.21
           85bp                         0.63                      0.63                          0.39                       0.21
           86bp                         0.58                      0.66                          0.33                       0.21
           87bp                         0.59                      0.63                          0.31                       0.20
           88bp                         0.55                      0.61                          0.30                       0.19
           89bp                         0.55                      0.62                          0.28                       0.19
           90bp                         0.58                      0.66                          0.33                       0.22
           91bp                         0.53                      0.62                          0.28                       0.20
           92bp                         0.58                      0.68                          0.31                       0.22
           93bp                         0.56                      0.63                          0.30                       0.20
           94bp                         0.58                      0.67                          0.32                       0.20
           95bp                         0.58                      0.62                          0.33                       0.19
           96bp                         0.57                      0.63                          0.31                       0.21
           97bp                         0.58                      0.62                          0.34                       0.19
           98bp                         0.52                      0.58                          0.27                       0.20
           99bp                         0.56                      0.59                          0.31                       0.19
          100bp                         0.53                      0.62                          0.27                       0.19
          120bp                         0.48                      0.59                          0.24                       0.18
          140bp                         0.52                      0.56                          0.25                       0.20
          160bp                         0.46                      0.36                          0.21                       0.20
          180bp                         0.48                      0.34                          0.23                       0.18
          200bp                         0.42                      0.34                          0.17                       0.18
          250bp                         0.41                      0.33                          0.15                       0.17
          300bp                         0.43                      0.33                          0.18                       0.18
          350bp                         0.40                      0.34                          0.15                       0.18
          400bp                         0.37                      0.32                          0.13                       0.18
          450bp                         0.39                      0.32                          0.14                       0.18
          500bp                         0.37                      0.30                          0.12                       0.18
          600bp                         0.36                      0.30                          0.11                       0.18
          700bp                         0.35                      0.30                          0.10                       0.18
          800bp                         0.34                      0.30                          0.09                       0.18
          900bp                         0.34                      0.29                          0.09                       0.19
         1000bp                         0.33                      0.29                          0.08                       0.18
         2000bp                         0.31                      0.27                          0.06                       0.18
         3000bp                         0.30                      0.27                          0.05                       0.18
         4000bp                         0.30                      0.27                          0.05                       0.18
         5000bp                         0.29                      0.27                          0.04                       0.18
         6000bp                         0.29                      0.27                          0.04                       0.18
         7000bp                         0.28                      0.26                          0.03                       0.18
"""
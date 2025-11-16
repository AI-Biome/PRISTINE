#!/usr/bin/env python3
"""
Similarity and distance metrics for sequence analysis.

This module contains various distance and similarity functions used for
comparing sequence profiles, k-mer frequencies, and other vector representations.
"""

from __future__ import annotations
import numpy as np
from collections import Counter
from datasketch import WeightedMinHashGenerator


def manhattan_distance(point1, point2):
    """
    Calculate the Manhattan distance between two points.

    :param point1: A list or tuple representing the coordinates of the first point.
    :param point2: A list or tuple representing the coordinates of the second point.
    :return: The Manhattan distance between the two points.
    """
    if len(point1) != len(point2):
        raise ValueError("Both points must have the same number of dimensions")

    distance = sum(abs(p1 - p2) for p1, p2 in zip(point1, point2))
    return distance


def cosine_similarity(vector_a, vector_b):
    """
    Calculates the cosine similarity between two vectors.

    :param vector_a: A list or tuple representing the first vector.
    :param vector_b: A list or tuple representing the second vector.
    :return: The cosine similarity between vector_a and vector_b.
    """
    # Ensure inputs are lists or tuples
    if not isinstance(vector_a, (list, tuple)) or not isinstance(vector_b, (list, tuple)):
        raise TypeError("Input vectors must be lists or tuples.")

    # Convert to NumPy arrays
    vector_a = np.array(vector_a)
    vector_b = np.array(vector_b)

    # Calculate dot product and magnitudes
    dot_product = np.dot(vector_a, vector_b)
    magnitude_a = np.linalg.norm(vector_a)
    magnitude_b = np.linalg.norm(vector_b)

    # Avoid division by zero
    if magnitude_a == 0 or magnitude_b == 0:
        return 0  # Return 0 for cosine similarity if either vector is zero

    # Compute cosine similarity
    return dot_product / (magnitude_a * magnitude_b)


def jaccard_similarity_with_frequencies(vector_a, vector_b):
    """
    Calculates the Jaccard similarity between two lists or tuples containing frequencies of individual items (e.g., p-mers).

    :param vector_a: A list or tuple representing the first vector.
    :param vector_b: A list or tuple representing the second vector.
    :return: The Jaccard similarity between the two collections, accounting for frequencies.
    """
    # Convert lists to Counters to handle frequencies
    counter_a = Counter(vector_a)
    counter_b = Counter(vector_b)

    # Calculate intersection and union based on minimum and maximum frequencies
    intersection = sum((min(counter_a[item], counter_b[item]) for item in counter_a if item in counter_b))
    union = sum((max(counter_a[item], counter_b[item]) for item in set(counter_a) | set(counter_b)))

    # Avoid division by zero in case both counters are empty
    if union == 0:
        return 0

    # Compute Jaccard similarity with frequencies
    return intersection / union


def calculate_weighted_minhash_similarity(frequency_vector_a, frequency_vector_b, num_perm=128):
    """
    Calculates the Weighted MinHash Jaccard similarity between two frequency vectors.

    :param frequency_vector_a: A list or tuple representing the frequency of p-mers for the first vector.
    :param frequency_vector_b: A list or tuple representing the frequency of p-mers for the second vector.
    :param num_perm: The number of permutations (hash functions) for the MinHash.
    :return: The Weighted MinHash Jaccard similarity between the two vectors.
    """
    if not isinstance(frequency_vector_a, (list, tuple)) or not isinstance(frequency_vector_b, (list, tuple)):
        raise TypeError("Input vectors must be lists or tuples.")

    if len(frequency_vector_a) != len(frequency_vector_b):
        raise ValueError("The two frequency vectors must be of the same length.")

    wmg = WeightedMinHashGenerator(len(frequency_vector_a), sample_size=num_perm)

    wm_a = wmg.minhash(frequency_vector_a)
    wm_b = wmg.minhash(frequency_vector_b)

    return wm_a.jaccard(wm_b)


def add_missing_keys_with_zero(dict1, dict2):
    """
    Adds keys that are present in dict1 but not in dict2 (and vice versa) to each dictionary with values set to zero.

    :param dict1: The first dictionary.
    :param dict2: The second dictionary.
    :return: Two dictionaries with added keys and values set to zero where keys were missing.
    """
    # Find keys unique to each dictionary
    keys_in_dict1_not_in_dict2 = dict1.keys() - dict2.keys()
    keys_in_dict2_not_in_dict1 = dict2.keys() - dict1.keys()

    # Add missing keys to each dictionary with values set to zero
    for key in keys_in_dict1_not_in_dict2:
        dict2[key] = 0
    for key in keys_in_dict2_not_in_dict1:
        dict1[key] = 0

    return dict1, dict2

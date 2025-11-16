#!/usr/bin/env python3
"""
Data structures for quasi-alignment and genomic segment representation.

This module contains classes for representing genomic segments and their
quasi-alignments, used in the QuasiAlignmentStrategy pipeline.
"""

from __future__ import annotations
import numpy as np
from similarity_metrics import manhattan_distance, add_missing_keys_with_zero


class Segment:
    """
    A class to represent an individual segment in a quasi-alignment. It includes information about the sequence ID, species,
    and the position of the segment within the sequence.
    """

    def __init__(self, seq_id: str, species_name: str, gene_name: str, start: int, length: int, sequence: str, p: int = 3):
        """
        Initialize a Segment object.

        :param seq_id: The ID of the sequence from which the segment comes.
        :param species_name: The species of the sequence.
        :param gene_name: The gene region of the sequence.
        :param start: The starting position of the segment within the sequence.
        :param length: The length of the segment within the sequence.
        :param sequence: The actual sequence string.
        :param p: The length of p-mers for composition analysis (default: 3).
        """
        self.seq_id = seq_id
        self.species_name = species_name
        self.gene_name = gene_name
        self.start = start      # 0-indexed
        self.length = length
        self.sequence = sequence
        self.pmer_profile = self.get_pmer_composition(self, p)

    def get_segment(self, sequence):
        """Returns the segment of the sequence based on start and length."""
        end = self.start + self.length
        return sequence[self.start:end]

    def get_pmer_composition(self, segment, p):
        """
        Takes a Segment object and computes the p-mer composition of the segment.
        The output is a dictionary where keys are only the p-mers present in the segment's sequence, and values are their counts.
        P-mers are stored in uppercase.

        :param segment: A Segment object containing the start and length information.
        :param p: The length of the p-mer (substring).
        :return: A dictionary with present p-mers as keys and their counts as values.
        """
        sequence = segment.sequence.upper()  # Convert the entire sequence to uppercase
        pmer_dict = {}

        # Loop through the sequence and extract p-mers
        for i in range(len(sequence) - p + 1):  # Ensure we don't go beyond the sequence
            pmer = sequence[i:i + p]

            # Add or update the count of the p-mer in the dictionary
            if pmer in pmer_dict:
                pmer_dict[pmer] += 1
            else:
                pmer_dict[pmer] = 1

        return pmer_dict

    def get_distance(self, other_segment, distance_func):
        """
        Calculates the distance between this segment and another segment using a specified distance function.

        :param other_segment: The other Segment object to calculate the distance to.
        :param distance_func: A function that takes two lists (or tuples) as input and returns the distance between them.
        :return: The distance between the two p-mer profiles.
        """
        if not callable(distance_func):
            raise ValueError("distance_func must be a callable function")

        aligned_profile_self, aligned_profile_other = add_missing_keys_with_zero(
            self.pmer_profile, other_segment.pmer_profile
        )

        profile_values_self = tuple(aligned_profile_self.values())
        profile_values_other = tuple(aligned_profile_other.values())

        # Calculate and return the distance using the specified distance function
        return distance_func(profile_values_self, profile_values_other)

    def __repr__(self):
        return f"Segment(seq_id={self.seq_id}, species={self.species_name}, start={self.start}, length={self.length})"


class QuasiAlignment:
    """
    A class to represent a quasi-alignment that contains multiple segments.
    Each segment is associated with a sequence from a particular species.
    """

    def __init__(self, cluster_id: str):
        """
        Initialize a QuasiAlignment object.

        :param cluster_id: A unique identifier for the quasi-alignment cluster.
        """
        self.cluster_id = cluster_id
        self.segments = []
        self.medoid = None

    def add_segment(self, new_segment: Segment):
        """
        Adds a new Segment to the quasi-alignment and recalculates the medoid.
        The medoid is the segment with the smallest average distance to all other segments.

        :param new_segment: Segment object to be added to the quasi-alignment.
        """
        # Add the new segment to the list
        self.segments.append(new_segment)

        # Recalculate the medoid if there is more than one segment
        if len(self.segments) > 1:
            # Calculate the average distance of each segment to all others
            min_avg_distance = float('inf')
            new_medoid = None

            for segment in self.segments:
                distances = [
                    segment.get_distance(other, manhattan_distance)
                    for other in self.segments if other != segment
                ]
                avg_distance = np.mean(distances)

                # Update the medoid if a smaller average distance is found
                if avg_distance < min_avg_distance:
                    min_avg_distance = avg_distance
                    new_medoid = segment

            # Set the new medoid
            self.medoid = new_medoid
        else:
            # If only one segment, it's the medoid by default
            self.medoid = new_segment

    def get_segments_by_species(self, species: str):
        """
        Get all segments that belong to a specific species.

        :param species: The species to filter segments by.
        :return: A list of segments that belong to the specified species.
        """
        return [segment for segment in self.segments if segment.species_name == species]

    def __repr__(self):
        return f"QuasiAlignment(cluster_id={self.cluster_id}, num_segments={len(self.segments)}, medoid={self.medoid})"

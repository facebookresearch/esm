from language import (
    FixedLengthSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
)


def free_hallucination(sequence_length: int) -> ProgramNode:
    sequence = FixedLengthSequenceSegment(sequence_length)
    return ProgramNode(
        sequence_segment=sequence,
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            MinimizeSurfaceHydrophobics(),
        ],
    )

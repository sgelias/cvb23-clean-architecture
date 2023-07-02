from typing import Callable, Generator, Literal

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from cvb23.core.domain.dtos.blast_config import BlastConfigDTO
from cvb23.core.domain.dtos.blast_results import BlastResultDataTable
from cvb23.core.domain.dtos.query_sequences import QuerySequencesDTO
from cvb23.core.domain.entities.blast_executor import BlastExecutor
from cvb23.core.domain.exceptions import SequentialBlastException
from cvb23.settings import ENV_NUM_THREADS

from .match_blast_type import match_blast_type


def run_sequential_blast(
    query_sequences: QuerySequencesDTO,
    blast_config: BlastConfigDTO,
    blast_executor: BlastExecutor,
) -> Literal[True]:
    """Run BLAST sequentially on a list of query sequences.

    Args:
        query_sequences (QuerySequencesDTO): Query sequences.
        blast_config (BlastConfigDTO): BLAST configuration.
        blast_executor (BlastExecutor): BLAST executor.

    Returns:
        list[BlastResultDTO]: List of BLAST results.

    Raises:
        SequentialBlastException: If an error occurs during BLAST execution.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? 1 - Validate entry parameters
        # ? --------------------------------------------------------------------

        if not isinstance(query_sequences, QuerySequencesDTO):
            raise TypeError("query_sequences must be of type QuerySequencesDTO")

        if not isinstance(blast_config, BlastConfigDTO):
            raise TypeError("blast_config must be of type BlastConfigDTO")

        if not isinstance(blast_executor, BlastExecutor):
            raise TypeError("blast_executor must be of type BlastExecutor")

        # ? --------------------------------------------------------------------
        # ? 2 - Try to match blast type
        # ? --------------------------------------------------------------------

        blast_type = match_blast_type(
            query_type=query_sequences.data_type,
            subject_type=blast_config.reference_database.data_type,
        )

        # ? --------------------------------------------------------------------
        # ? 3 - Execute BLAST
        # ? --------------------------------------------------------------------

        blast_chunk_size = len(query_sequences) / ENV_NUM_THREADS

        chunk_function, sequences = __chunk_query_sequences(
            query=query_sequences,
            chunk_size=blast_chunk_size,
        )

        blast_results_list: list[BlastResultDataTable] = []
        for chunk in chunk_function(sequences):
            blast_results_list.extend(
                blast_executor.run(
                    query_sequences=chunk,
                    blast_config=blast_config,
                    blast_type=blast_type,
                )
            )

        # ? --------------------------------------------------------------------
        # ? 4 - Evaluate BLAST results
        # ? --------------------------------------------------------------------

        concatenated_blast_results: pd.DataFrame = (
            BlastResultDataTable.validate(pd.concat(blast_results_list, axis=1))
        )

        grouped_blast_results = (
            concatenated_blast_results.groupby([BlastResultDataTable.query])
            .apply(
                lambda x: x.sort_values(
                    by=[BlastResultDataTable.bit_score],
                    ascending=False,
                )
            )
            .max()
        )

        grouped_blast_results.to_csv(
            query_sequences.path.parent.joinpath("blast").with_suffix(".tsv"),
        )

        # ? --------------------------------------------------------------------
        # ? 5 - Return a positive response
        # ? --------------------------------------------------------------------

        return True

    except Exception as e:
        raise SequentialBlastException(e)


def __chunk_query_sequences(
    query: QuerySequencesDTO,
    chunk_size: int,
) -> tuple[
    Callable[[list[SeqRecord]], Generator[list[SeqRecord], None, None]],
    list[SeqRecord],
]:
    """Chunk query sequences.

    Args:
        query (QuerySequencesDTO): Query sequences.
        chunk_size (int): Chunk size.

    Returns:
        tuple[
            Callable[[list[SeqRecord]], Generator[list[SeqRecord], None, None]],
            list[SeqRecord]
        ]:
            Chunked query sequences.

    """

    def __chunks(
        query_sequences: list[SeqRecord],
    ) -> Generator[list[SeqRecord], None, None]:
        """Chunk sequences list generator.

        Args:
            query_sequences (list[SeqRecord]): Query sequences.

        Yields:
            Generator[list[SeqRecord], None, None]: Chunked query sequences.

        """

        for i in range(0, len(query_sequences), chunk_size):
            yield query_sequences[i : i + chunk_size]

    try:
        return (__chunks, [i for i in SeqIO.parse(query.path, query.format)])

    except Exception as e:
        raise SequentialBlastException(e)
